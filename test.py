from equilibrator_api import ComponentContribution, Reaction
import numpy as np
import pandas as pd
import cobra
from cobra.core.reaction import Reaction as cobraReaction
from cobra.flux_analysis.variability import flux_variability_analysis
from cobra.util.solver import set_objective
import cvxopt
from cvxopt import glpk
import networkx as nx
import json
import time
import os
from six import iteritems
import re

# Constants (dG0 in kJ/mol) original dG0 data units
R = 8.3144598 * 1e-6 # kJ. K^-1. mmol^-1
T = 310.15 # 298.15 # K

# Parameters
loopless_fva = False
numerator_metabolites = p_mets = None
denominator_metabolites = q_mets = None
uncertainty_threshold = alpha = 100
biomass_threshold = beta = 0.5
dG0_error_fraction = gamma = 1
Gibbs_eps = 1e-8 # kJ/mmol
M = 1e8
x_min, x_max = 1e-4, 2e2 # mM
Intracellular_pH = pH_i = 7.3
pH_deviation = delta_pH = 0.3
# Fluxes in mmol.gDW^-1.min^-1
# (in Teppe et al 2013* they use 1e-5, 1e-1!
# https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0075370)
fileName = 'graphData.json'
dirName = 'Cond2'
notes = 'Growth on Glucose, 10 mmol.min^-1.gDW^-1, MILP formulation!'

def updateExchangeReactionBounds(GEM, uptakeRxn='EX_glc__D_e_backward', carbonUptakeRate=10):
    """
    Update exchange reaction bounds to simulate appropriate growth medium
    conditions.
    """
    def isOrganicExchange(ID):
        compoundAtoms = list(GEM.reactions.get_by_id(ID).reactants[0].formula)
        cond = (('C' in compoundAtoms)
                & ('H' in compoundAtoms)
                & ('o' not in compoundAtoms))  # discards cobalamine
        return cond

    ExchangeRxnIDs = [rxn.id for rxn in GEM.exchanges if len(rxn.products) == 0]
    for ID in ExchangeRxnIDs:
        try:
            if isOrganicExchange(ID):
                GEM.reactions.get_by_id(ID).lower_bound = 0
        except Exception:
            pass
    GEM.reactions.get_by_id(uptakeRxn).lower_bound = -carbonUptakeRate


def convert_to_irreversible_SRE(cobra_model):
    """Split reversible reactions into two irreversible reactions

    These two reactions will proceed in opposite directions. This
    guarentees that all reactions in the model will only allow
    positive flux values, which is useful for some modeling problems.

    cobra_model: A Model object which will be modified in place.
    Modified by Semidan Robaina, February 2019.
    Only splits non-blocked reactions under the bounds specified in the model

    """
    reactions_to_add = []
    coefficients = {}
    def onlyBackward(reaction):
        return reaction.lower_bound < 0 and reaction.upper_bound <= 0

    def backwardAndForward(reaction):
        return reaction.lower_bound < 0 and reaction.upper_bound > 0

    def changeDirection(reaction):
        def swapSign(number):
            return -number
        lb = swapSign(reaction.upper_bound)
        ub = swapSign(reaction.lower_bound)
        reaction.lower_bound = lb
        reaction.upper_bound = ub
        reaction.objective_coefficient * -1
        reaction.notes["reflection"] = 'reversed, only backward'
        for met in reaction._metabolites.keys():
            reaction._metabolites[met] *= -1

    def createBackwardReaction(reaction):
        backward_reaction = cobraReaction(reaction.id + '_backward')
        backward_reaction.lower_bound = 0
        backward_reaction.upper_bound = -reaction.lower_bound
        reaction_dict = {k: v * -1
                         for k, v in iteritems(reaction._metabolites)}
        backward_reaction.add_metabolites(reaction_dict)
        backward_reaction._model = reaction._model
        backward_reaction._genes = reaction._genes
        for gene in reaction._genes:
            gene._reaction.add(backward_reaction)
        backward_reaction.subsystem = reaction.subsystem
        backward_reaction.name = reaction.name + '_backward'
        backward_reaction._gene_reaction_rule = reaction._gene_reaction_rule
        coefficients[backward_reaction] = reaction.objective_coefficient * -1
        return backward_reaction

    for reaction in cobra_model.reactions:
        if onlyBackward(reaction):
            changeDirection(reaction)
        elif backwardAndForward(reaction):
            backward_reaction = createBackwardReaction(reaction)
            reactions_to_add.append(backward_reaction)
            reaction.id += '_forward'
            reaction.name += '_forward'
            reaction.lower_bound = 0

    cobra_model.add_reactions(reactions_to_add)
    set_objective(cobra_model, coefficients, additive=True)


# Load model
GEM = cobra.io.load_json_model('iJO1366.json')
updateExchangeReactionBounds(GEM, uptakeRxn='EX_glc__D_e', carbonUptakeRate=10)

# Find non-default irreversible reactions
print('Running flux variability analysis...')
fva = flux_variability_analysis(GEM, reaction_list=None,
                                fraction_of_optimum=beta, processes=2,
                                loopless=loopless_fva).round(decimals=8)
# fva.to_csv('iJO1366fva.csv')
GEM.optimize()
biomass_reaction = GEM.reactions.get_by_id('BIOMASS_Ec_iJO1366_core_53p95M')
biomass_reaction.lower_bound = beta * GEM.objective.value

# Update lower bound of fixed reactions and eliminate blocked ones
for rxn_id in fva.index:
    v_min, v_max = fva.loc[rxn_id].minimum, fva.loc[rxn_id].maximum
    if v_min > 0 or v_max < 0:
        GEM.reactions.get_by_id(rxn_id).lower_bound = v_min
        GEM.reactions.get_by_id(rxn_id).upper_bound = v_max
    if v_min == 0 and v_max == 0:
        GEM.reactions.get_by_id(rxn_id).remove_from_model(remove_orphans=True)

convert_to_irreversible_SRE(GEM)
GEM.reactions.get_by_id('EX_glc__D_e')

# Load equilibrator data
eq_api_7 = ComponentContribution(pH=7, ionic_strength=0.25)
eq_api_75 = ComponentContribution(pH=7.5, ionic_strength=0.25)

# Get iJO1366 reactions as KEGG reaction strings
GEM2KEGG = pd.read_csv(
    'equilibrator-api/src/equilibrator_api/data/iJO1366_reactions.csv',
                           index_col=0)

rxn = Reaction.parse_formula(GEM2KEGG.iloc[20].item())
dG0_prime_7, dG0_uncertainty_7 = np.array(eq_api_7.dG0_prime(rxn))
dG0_prime_75, dG0_uncertainty_75 = np.array(eq_api_75.dG0_prime(rxn))

print(dG0_prime_7, dG0_uncertainty_7, dG0_prime_75, dG0_uncertainty_75)

pH_i = 7.3
delta_pH = 0.2
np.log(10**(-(pH_i + delta_pH) + 3))
10**-6
10**(-(pH_i - delta_pH) + 3)
