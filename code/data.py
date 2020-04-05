import numpy as np
import pandas as pd
from six import iteritems
import re
import cobra
from cobra.core.reaction import Reaction as cobraReaction
from cobra.flux_analysis.variability import flux_variability_analysis
from cobra.util.solver import set_objective
from equilibrator_api import ComponentContribution, Q_, parse_reaction_formula


def updateExchangeReactionBounds(GEM,
                                 uptakeRxn='EX_glc__D_e',
                                 carbonUptakeRate=10):
    """
    Update exchange reaction bounds to simulate appropriate growth medium
    conditions. Default max glucose uptake rate in iML1515 is 10 mmol/(min.gDW)
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
        reaction.objective_coefficient *= -1
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
            reaction.lower_bound = 0
            reaction.id += '_forward'
            reaction.name += '_forward'

    GEM.add_reactions(reactions_to_add)
    set_objective(GEM, coefficients, additive=True)


def removeBlockedReactions(GEM, v_eps=None):        
    blockedRxns = cobra.flux_analysis.find_blocked_reactions(
        GEM, zero_cutoff=v_eps, open_exchanges=False)
    for rxn in blockedRxns:
        GEM.reactions.get_by_id(rxn).remove_from_model(remove_orphans=True)
        
        
def removeBlockedReactionsLoopless(GEM):
    fva = flux_variability_analysis(GEM, reaction_list=None,
                                    fraction_of_optimum=0, processes=2,
                                    loopless=True).round(decimals=8)
    for rxn_id in fva.index:
        v_min, v_max = fva.loc[rxn_id].minimum, fva.loc[rxn_id].maximum
        if v_min == 0 and v_max == 0:
            GEM.reactions.get_by_id(rxn_id).remove_from_model(remove_orphans=True)
            
            
def findReactionsWithPositiveFluxAtOptimum(GEM, fraction_of_optimum=1):
    """from cobra.flux_analysis.variability import flux_variability_analysis
    GEM must be irreversible, i.e., split reversible reactions
    """
    positive_rxns = []
    fva = flux_variability_analysis(GEM, reaction_list=None,
                                    fraction_of_optimum=fraction_of_optimum, 
                                    processes=2, loopless=False).round(decimals=8)
    for rxn_id in fva.index:
        v_min, v_max = fva.loc[rxn_id].minimum, fva.loc[rxn_id].maximum
        if v_min > 0:
            positive_rxns.append(rxn_id)
    return positive_rxns


def findMetaboliteByCommonName(GEM, met_name):
    candidates = []
    met_ids = [met.id for met in GEM.metabolites]
    met_names = [met.name for met in GEM.metabolites]
    for i, met_id in enumerate(met_ids):
        if (met_name.lower() in met_id.lower() or
            met_name.lower() in met_names[i].lower()):
            candidates.append([met_id, met_names[i]])
    return candidates


def findMetabolitesInPathways(GEM, systems_df, pathways=['Carbohydrate metabolism']):
    """
    Find metabolites participating in reactions of central carbon metabolism
    """
    systems = systems_df['Subsystems'].loc[
        systems_df['Systems'].isin(pathways)].values
    
    def isPathwayMetabolite(met):
        met_subsystems = [rxn.subsystem for rxn in met.reactions]
        return len(np.intersect1d(systems, met_subsystems)) > 0
        
    pathway_metabolites = []
    for met in GEM.metabolites:
        if isPathwayMetabolite(met):
            pathway_metabolites.append(met.id)
    return pathway_metabolites


def findThermodynamicallyInfeasibleReactions(GEM, dG0data, x_min, x_max, 
                                             fraction_of_optimum=1):
    """
    Find reactions that are needed to carry flux to produce biomass but that
    are thermodynamically infeasible due to concentration and dG0 constraints
    """
    infeasible_rxns = {}
    positive_rxns = findReactionsWithPositiveFluxAtOptimum(GEM, 
                                                           fraction_of_optimum=fraction_of_optimum)
    
    dG0_lb = {}
    for rxn_id in dG0data.keys():
        dG0_lb[rxn_id] = dG0data[rxn_id]['dG0'] - dG0data[rxn_id]['error']
        
        
    for rxn_id in dG0data.keys():
        dG0_lb = dG0data[rxn_id]['dG0'] - dG0data[rxn_id]['error']
        rxn = GEM.reactions.get_by_id(rxn_id)
        
        logx_sum = 0
        for met in rxn.metabolites.keys():
            coeff = rxn.metabolites[met]
            extreme_x = 0
            if coeff < 0:
                extreme_x = x_max[met.id]
            else:
                extreme_x = x_min[met.id]
            logx_sum += coeff * np.log(extreme_x)
        
        dG0_r = dG0_lb + par.R * par.T * logx_sum
        if dG0_r > 0:
            infeasible_rxns[rxn_id] = {'dG0_r':dG0_r,'dG0_lb': dG0_lb,
                                       'sum_logx':par.R * par.T * logx_sum}
            
    return infeasible_rxns
        

def prepareGEM(path_to_GEM, carbon_source=None,
               uptake_rate=None, loopless=False,
               biomass_reaction_id=None):

    GEM = cobra.io.load_json_model(path_to_GEM)
    updateExchangeReactionBounds(GEM, uptakeRxn=carbon_source,
                                 carbonUptakeRate=uptake_rate)
    
    if loopless:
        removeBlockedReactionsLoopless(GEM)
    else:
        removeBlockedReactions(GEM, v_eps=None)
    
    convert_to_irreversible_SRE(GEM)
    _ = cobra.manipulation.delete.prune_unused_metabolites(GEM)

    return GEM


def getFreeEnergyData(GEM, work_directory='', pH_i=None,
                      Ionic_strength=None,
                      dG0_uncertainty_threshold=None):

    # Load equilibrator data
    eq_api = ComponentContribution(p_h=Q_(str(pH_i)), ionic_strength=Q_(Ionic_strength))

    # Get iJO1366/iML1515 reactions as KEGG reaction strings
    GEM2KEGG = pd.read_csv(work_directory + '/iJO1366_reactions.csv', index_col=0)

    # Prepare dictionary with dG0
    def originallyIrreversible(reaction_id):
        return 'forward' not in reaction_id and 'backward' not in reaction_id

    dG0_data = {}
    GEM_rxns = np.array([rxn.id for rxn in GEM.reactions])

    for rxn_id in GEM_rxns:
        try:
            if originallyIrreversible(rxn_id):
                id, direction = rxn_id, ''
            else:
                id, direction = re.split('_(forward|backward)', rxn_id)[:2]

            rxn = parse_reaction_formula(GEM2KEGG.loc[id.lower()]['formula'])

            dG0_prime = eq_api.standard_dg_prime(rxn).value.magnitude  # kJ/mmol
            dG0_uncertainty = eq_api.standard_dg_prime(rxn).error.magnitude  # kJ/mol
            if dG0_uncertainty < abs(dG0_uncertainty_threshold * dG0_prime):
                # remove uncertain dG0 data
                if 'backward' in direction:
                    dG0_data[rxn_id] = {'dG0': -dG0_prime,
                                        'error': dG0_uncertainty}
                else:
                    dG0_data[rxn_id] = {'dG0': dG0_prime,
                                        'error': dG0_uncertainty}
        except Exception:
            pass

    return dG0_data
