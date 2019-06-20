import numpy as np
import pandas as pd
import json
import time
import os
from six import iteritems
import re

import cobra
from cobra.core.reaction import Reaction as cobraReaction
from cobra.flux_analysis.variability import flux_variability_analysis
from cobra.util.solver import set_objective
from equilibrator_api import ComponentContribution, Reaction
import cvxopt
from cvxopt import glpk
import networkx as nx


# from parameters import * #that would solve the problem
def updateExchangeReactionBounds(GEM,
                                 uptakeRxn='EX_glc__D_e',
                                 carbonUptakeRate=10):
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
    # Temporarily remove constraint on ATPM
    GEM.reactions.get_by_id('ATPM').lower_bound = 0
    return GEM

def convert_to_irreversible_SRE(GEM):
    """Split reversible reactions into two irreversible reactions: one going in
    the forward direction, the other in the backward direction. In this manner,
    all reactions in the model carry non-negative flux values. Forward reactions
    are tagged as "forward" while backward reactions as "backward".

    GEM: A Model object which will be modified in place.
    Modified from the deprecated cobrapy version by Semidan Robaina,
    February 2019.
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

    for reaction in GEM.reactions:
        if onlyBackward(reaction):
            changeDirection(reaction)
        elif backwardAndForward(reaction):
            backward_reaction = createBackwardReaction(reaction)
            reactions_to_add.append(backward_reaction)
            reaction.id += '_forward'
            reaction.name += '_forward'
            reaction.lower_bound = 0

    GEM.add_reactions(reactions_to_add)
    set_objective(GEM, coefficients, additive=True)

    return GEM


def prepareGEM(work_directory='', carbon_source=None,
                uptake_rate=None, biomass_threshold=None,
                loopless_fva=None,
                biomass_reaction_id=None):

    # Load model (Update exchange reaction bounds)
    GEM = cobra.io.load_json_model(work_directory + 'iJO1366.json')
    updateExchangeReactionBounds(GEM, uptakeRxn=carbon_source, carbonUptakeRate=uptake_rate)

    # Find non-default irreversible reactions
    print('Running flux variability analysis...')
    fva = flux_variability_analysis(GEM, reaction_list=None,
                                    fraction_of_optimum=biomass_threshold,
                                    processes=2,
                                    loopless=loopless_fva).round(decimals=8)
    # fva.to_csv('iJO1366fva.csv')
    GEM.optimize()
    biomass_reaction = GEM.reactions.get_by_id(biomass_reaction_id)
    biomass_reaction.lower_bound = biomass_threshold * GEM.objective.value

    # Update lower bound of fixed reactions and eliminate blocked ones
    for rxn_id in fva.index:
        v_min, v_max = fva.loc[rxn_id].minimum, fva.loc[rxn_id].maximum
        if v_min > 0 or v_max < 0:
            GEM.reactions.get_by_id(rxn_id).lower_bound = v_min
            GEM.reactions.get_by_id(rxn_id).upper_bound = v_max
        if v_min == 0 and v_max == 0:
            GEM.reactions.get_by_id(rxn_id).remove_from_model(remove_orphans=True)

    # Convert model to irreversible
    convert_to_irreversible_SRE(GEM)

    # Remove unused metabolites
    _ = cobra.manipulation.delete.prune_unused_metabolites(GEM)

    return GEM


def getFreeEnergyData(GEM, work_directory='', pH_i=None,
                      Ionic_strength=None,
                      dG0_uncertainty_threshold=None):

    # Load equilibrator data
    eq_api = ComponentContribution(pH=pH_i, ionic_strength=Ionic_strength)

    # Get iJO1366 reactions as KEGG reaction strings
    GEM2KEGG = pd.read_csv(work_directory
         + 'equilibrator-api/src/equilibrator_api/data/iJO1366_reactions.csv',
            index_col=0)

    # Prepare dictionary with dG0
    def originallyIrreversible(reaction_id):
        return 'forward' not in reaction_id and 'backward' not in reaction_id

    dG0_data = {}
    GEM_rxns = np.array([rxn.id for rxn in GEM.reactions])
    GEM_mets = np.array([met.id for met in GEM.metabolites])

    for rxn_id in GEM_rxns:
        try:
            if originallyIrreversible(rxn_id):
                id, direction = rxn_id, ''
            else:
                id, direction = re.split('_(forward|backward)', rxn_id)[:2]

            rxn = Reaction.parse_formula(GEM2KEGG.loc[id.lower()].item())
            dG0_prime, dG0_uncertainty = np.array(eq_api.dG0_prime(rxn)) * 1e-3 # convert to kJ/mmol
            if dG0_uncertainty < abs(dG0_uncertainty_threshold * dG0_prime):
                # remove uncertain dG0 data
                if 'backward' in direction:
                    dG0_data[rxn_id] = [-dG0_prime, dG0_uncertainty]
                else:
                    dG0_data[rxn_id] = [dG0_prime, dG0_uncertainty]
        except Exception:
            pass

    return dG0_data
