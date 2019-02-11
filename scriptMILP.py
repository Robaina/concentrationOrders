if __name__ == '__main__':

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

# Constants (dG0 in kJ/mol)
R = 8.3144598 * 1e-3 # kJ. K^-1. mol^-1
T = 310.15 # 298.15 # K

# Parameters
numerator_metabolites = p_mets = ['h_p']
denominator_metabolites = q_mets = ['h_c']
uncertainty_threshold = alpha = 0.5
biomass_threshold = beta = 1
dG0_error_fraction = gamma = 1
Gibbs_eps = 1e-6 # kJ/mol
M = 1e6
x_min, x_max = 1e-6, 1e-1
# (in Teppe et al 2013* they use 1e-5, 1e-1!
# https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0075370)
fileName = 'graphData.json'
dirName = 'Cond1'
notes = "Original iJO1366 flux bounds, loopless=True, remove_orphans=True, in nanomolar"

# ******************************* Load model, dG0********************************
# *******************************************************************************

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

    def changeReactionDirection(reaction):
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
        backward_reaction = cobraReaction(reaction.id + '_backward') #check correct path
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
            changeReactionDirection(reaction)
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

# Find non-default irreversible reactions
print('Running flux variability analysis...')
fva = flux_variability_analysis(GEM, reaction_list=None,
                                fraction_of_optimum=beta, processes=2,
                                loopless=False).round(decimals=4)
fva.to_csv('iJO1366fva.csv')
GEM.reactions[6].lower_bound = fva['minimum'].loc['BIOMASS_Ec_iJO1366_core_53p95M']

for rxn_id in fva.index:
    v_min, v_max = fva.loc[rxn_id].minimum, fva.loc[rxn_id].maximum
    if v_min < 0 and v_max <= 0:
        ub = min(0, GEM.reactions.get_by_id(rxn_id).upper_bound)
        GEM.reactions.get_by_id(rxn_id).upper_bound = ub
    elif v_min >= 0 and v_max > 0:
        lb = max(0, GEM.reactions.get_by_id(rxn_id).lower_bound)
        GEM.reactions.get_by_id(rxn_id).lower_bound = lb
# fva['minimum'].loc['EX_glc__D_e']
# GEM_rxns
# GEM.reactions.get_by_id('EX_glc__D_e')

# Remove blocked reactions (using fva bounds)
blockedRxns = cobra.flux_analysis.find_blocked_reactions(
    GEM, zero_cutoff=1e-9, open_exchanges=False)
for rxn in blockedRxns:
    GEM.reactions.get_by_id(rxn).remove_from_model(remove_orphans=True)

convert_to_irreversible_SRE(GEM)
_ = cobra.manipulation.delete.prune_unused_metabolites(GEM)

# Load equilibrator data
eq_api = ComponentContribution(pH=7.5, ionic_strength=0.25)

# Get iJO1366 reactions as KEGG reaction strings
GEM2KEGG = pd.read_csv(
    'equilibrator-api/src/equilibrator_api/data/iJO1366_reactions.csv',
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
        dG0_prime, dG0_uncertainty = eq_api.dG0_prime(rxn)
        if 'backward' in direction:
            dG0_data[rxn_id] = [-dG0_prime, dG0_uncertainty]
        else:
            dG0_data[rxn_id] = [dG0_prime, dG0_uncertainty]
    except Exception:
        pass

# ************************************Prepare MILP************************************
# ************************************************************************************

# Categorize reactions
def convertRxnIDtoIndex(GEM, rxn_id):
    return [rxn.id for rxn in GEM.reactions].index(rxn_id)
def convertMetIDtoIndex(GEM, met_id):
    return [met.id for met in GEM.metabolites].index(met_id)
def convertRxnIndexToID(GEM, rxn_index):
    return [rxn.id for rxn in GEM.reactions][rxn_index]
def convertMetIndexToID(GEM, met_index):
    return [met.id for met in GEM.metabolites][met_index]

Irr_rxns_with_dG0 = [rxn.id for rxn in GEM.reactions
                     if ('forward' not in rxn.id and 'backward' not in rxn.id
                         and rxn.id in dG0_data.keys())]
For_rxns_with_dG0 = [rxn.id for rxn in GEM.reactions
                     if ('forward' in rxn.id and rxn.id in dG0_data.keys())]
Back_rxns_with_dG0 = [rxn.id for rxn in GEM.reactions
                     if ('backward' in rxn.id and rxn.id in dG0_data.keys())]
Rest_rxns = [rxn.id for rxn in GEM.reactions
             if (rxn.id not in Irr_rxns_with_dG0
                 and rxn.id not in For_rxns_with_dG0
                 and rxn.id not in Back_rxns_with_dG0)]

Irr_rxns_with_dG0_Idx = [convertRxnIDtoIndex(GEM, id) for id in Irr_rxns_with_dG0]
For_rxns_with_dG0_Idx = [convertRxnIDtoIndex(GEM, id) for id in For_rxns_with_dG0]
Back_rxns_with_dG0_Idx = [convertRxnIDtoIndex(GEM, id) for id in Back_rxns_with_dG0]
Rest_rxns_Idx = [convertRxnIDtoIndex(GEM, id) for id in Rest_rxns]

# Extract stoichiometric matrix
N_original_order = cobra.util.array.create_stoichiometric_matrix(GEM, array_type='dense')
N_mets, N_rxns = N_original_order.shape

# Re-order reactions: Irr_rxns_with_dG0, For_rxns_with_dG0, Back_rxns_with_dG0, Rest_rxns
re_ordered_rxns = np.concatenate((Irr_rxns_with_dG0, For_rxns_with_dG0,
                                     Back_rxns_with_dG0, Rest_rxns))
re_ordered_indices = [convertRxnIDtoIndex(GEM, id) for id in re_ordered_rxns]
rxns_with_dG0 = np.concatenate((Irr_rxns_with_dG0, For_rxns_with_dG0,
                                     Back_rxns_with_dG0))
rxns_with_dG0_indices = [convertRxnIDtoIndex(GEM, id) for id in rxns_with_dG0]

N = N_original_order[:, re_ordered_indices]
N_met_pairs = int(0.5 * N_mets * (N_mets - 1))

print(('There are ' + str(N_rxns - len(Rest_rxns))
       + ' reactions with dG0 data and '
       + str(N_met_pairs) + ' metabolite pairs'))

# Standard in cvxopt is Ax <= b so have to change signs and add epsilon to rhs
# Bounds
logx_min = (np.log(x_min)) * np.ones((N_mets, 1))
logx_max = (np.log(x_max)) * np.ones((N_mets, 1))
v_min = np.array([rxn.lower_bound for rxn in GEM.reactions])

v_max = np.array([rxn.upper_bound for rxn in GEM.reactions])
# v_min.shape, v_max.shape = (N_rxns, 1), (N_rxns, 1)

# Construct vectors: dG0min, dG0max
N_rxns_dG0 = len(rxns_with_dG0)
dG0min, dG0max = np.zeros((N_rxns_dG0, 1)), np.zeros((N_rxns_dG0, 1))
for i, rxn in enumerate(rxns_with_dG0):
    dG0_i, dG0_error_i = dG0_data[rxn]
    dG0min[i] = dG0_i - gamma * dG0_error_i
    dG0max[i] = dG0_i + gamma * dG0_error_i

# Define dimensions of submatrices
N_rev_with_dG0 = len(For_rxns_with_dG0)
N_rxns_no_dG0 = len(Rest_rxns)
N_irr_with_dG0 = len(Irr_rxns_with_dG0)

# Build constraints matrix (cols: N_mets + 2*N_rxns_dG0 + N_rxns )
A0 = np.hstack((np.zeros((N_mets, N_mets + 2*N_rxns_dG0)), N)) #Nv = 0

A1 = np.hstack((R*T*N[:, rxns_with_dG0_indices].transpose(), np.identity(N_rxns_dG0),
                M * np.identity(N_rxns_dG0), np.zeros((N_rxns_dG0, N_rxns))))

A2 = np.hstack((np.zeros((N_rxns_dG0, N_mets)), -np.identity(N_rxns_dG0),
                np.zeros((N_rxns_dG0, N_rxns_dG0 + N_rxns))))

A3 = np.hstack((np.zeros((N_rxns_dG0, N_mets)), np.identity(N_rxns_dG0),
                np.zeros((N_rxns_dG0, N_rxns_dG0 + N_rxns))))

A4 = np.hstack((np.zeros((N_rev_with_dG0, N_mets + N_rxns_dG0 + N_irr_with_dG0)),
                np.diag(v_min[For_rxns_with_dG0_Idx]),
                np.zeros((N_rev_with_dG0, N_rev_with_dG0 + N_irr_with_dG0)),
                -np.identity(N_rev_with_dG0), np.zeros((N_rev_with_dG0,
                                                        N_rev_with_dG0 + N_rxns_no_dG0))))

A5 = np.hstack((np.zeros((N_rev_with_dG0, N_mets + N_rxns_dG0 + N_irr_with_dG0 + N_rev_with_dG0)),
                np.diag(-v_max[For_rxns_with_dG0_Idx]),
                np.zeros((N_rev_with_dG0, N_irr_with_dG0 + N_rev_with_dG0)),
                np.identity(N_rev_with_dG0), np.zeros((N_rev_with_dG0, N_rxns_no_dG0))))

A6 = np.hstack((np.zeros((N_rev_with_dG0, N_mets + N_rxns_dG0 + N_irr_with_dG0)),
                np.identity(N_rev_with_dG0), np.identity(N_rev_with_dG0),
                np.zeros((N_rev_with_dG0, N_rxns))))

A7 = np.hstack((-np.identity(N_mets), np.zeros((N_mets, 2*N_rxns_dG0 + N_rxns)))) #logx

A8 = np.hstack((np.identity(N_mets), np.zeros((N_mets, 2*N_rxns_dG0 + N_rxns)))) #logx

A9 = np.hstack((np.zeros((N_rxns_no_dG0, N_mets + 3*N_rxns_dG0)),
               -np.identity(N_rxns_no_dG0))) #v no dG0

A10 = np.hstack((np.zeros((N_rxns_no_dG0, N_mets + 3*N_rxns_dG0)),
                np.identity(N_rxns_no_dG0))) #v no dG0

A = cvxopt.matrix(np.vstack((A0, A1, A2, A3, A4, A5, A6, A7, A8, A9, A10)))

b = cvxopt.matrix(np.vstack((np.zeros((N_mets, 1)),
                             -Gibbs_eps * np.ones((N_rxns_dG0, 1)) + M,
                             -dG0min, dG0max,
                             np.zeros((2 * N_rev_with_dG0, 1)),
                             2 * np.ones((N_rev_with_dG0, 1)),
                             -logx_min, logx_max,
                             -v_min[Rest_rxns_Idx].reshape((-1, 1)),
                              v_max[Rest_rxns_Idx].reshape((-1, 1)))))

BinaryVariables = set(range(N_mets + N_rxns_dG0 + 1, N_mets + 2 * N_rxns_dG0 + 1))

# *************************Loop over metabolites******************************************
# ****************************************************************************************

# Start iteration over metabolite pairs (variable components)
met_orders = []
if p_mets is None:
    p_mets = range(N_mets)
if q_mets is None:
    q_mets = range(N_mets)
p_mets = ['fadh2_c']
q_mets = ['fad_c']
if p_mets is not None:
    p_mets = [convertMetIDtoIndex(GEM, met_id) for met_id in p_mets]
if q_mets is not None:
    q_mets = [convertMetIDtoIndex(GEM, met_id) for met_id in q_mets]

for p in p_mets:
    for q in q_mets:
        if p != q:
            # Find minimum ratio
            c = np.zeros(N_mets + 2 * N_rxns_dG0 + N_rxns)
            c[[p, q]] = [1, -1]
            c = cvxopt.matrix(c)

            res = glpk.ilp(c, A, b, B=set(BinaryVariables),
                                    options={'glpk':{'msg_lev':'GLP_MSG_OFF'}})
            # res = cvxopt.solvers.lp(c, A, b, solver='glpk',
            #                         options={'glpk':{'msg_lev':'GLP_MSG_OFF'}})
            x = res[1]
            z = x[p] - x[q]
            if z > 0:
                met_i, met_j = GEM_mets[[p, q]]
                min_ratio = np.e**z

                # Find maximum ratio
                c *= -1
                res = glpk.ilp(c, A, b, B=set(BinaryVariables),
                                        options={'glpk':{'msg_lev':'GLP_MSG_OFF'}})

                x = res[1]
                z = -x[p] + x[q]
                max_ratio = np.e**z

                met_orders.append([met_i, met_j, min_ratio, max_ratio])

z
met_orders
res
v_min[5]
x = res[1]
[x[i] for i in range(N_mets)]
x[p_mets[0]] - x[q_mets[0]]
def getMILPIndexFromReactionID(rxn_id):
    return N_mets + 2*N_rxns_dG0 + re_ordered_rxns.tolist().index(rxn_id)
x[getMILPIndexFromReactionID('BIOMASS_Ec_iJO1366_WT_53p95M')]
getMILPIndexFromReactionID('BIOMASS_Ec_iJO1366_WT_53p95M')
x[N_mets + 2*N_rxns_dG0 + 1]
min(v_max)

getMILPIndexFromReactionID('BIOMASS_Ec_iJO1366_WT_53p95M')

GEM.reactions.get_by_id('BIOMASS_Ec_iJO1366_WT_53p95M')

convertRxnIndexToID(GEM, 5)
[rxn.id for rxn in GEM.reactions]
GEM.reaction
res
from cvxopt import glpk
glpk

# SOME PROBLEM with the ILP solver, check it out!

    # ************************************* Graph************************************
    # *******************************************************************************

    # Build full graph adjacency matrix
    total_mets = np.unique(
        np.array(
            [[met_pairs[0], met_pairs[1]]
             for met_pairs in met_orders]).flatten()).tolist()
    N_nodes = len(total_mets)

    A, B = np.zeros((N_nodes, N_nodes)), np.zeros((N_nodes, N_nodes))
    for met_pair in met_orders:
        i, j = total_mets.index(
            met_pair[0]), total_mets.index(met_pair[1])
        A[i, j] = met_pair[2]
        B[i, j] = met_pair[3]

    # Get transitive reduction of the graph
    G = nx.from_numpy_matrix(A, create_using=nx.DiGraph())
    G_plus = nx.transitive_reduction(G)

    # Get adjacency matrix of the transitive reduction
    A_plus = nx.to_numpy_matrix(G)

    # Get list of edges
    reduced_met_orders = []
    for i in range(N_nodes):
        for j in range(N_nodes):
            min_ratio = A_plus[i, j]
            if min_ratio > 1:
                max_ratio = B[i, j]
                reduced_met_orders.append(
                     [total_mets[i], total_mets[j], min_ratio, max_ratio])

    # Build data dictionary
    data = {}
    data['nodes'] = []
    data['edges'] = []

    for met in total_mets:
        data['nodes'].append(
        {
            'data': {
                'id': met,
                'label': met,
                'name': iJO1366.metabolites.get_by_id(met).name
            }
        })

    for met_pairs in reduced_met_orders:
        data['edges'].append(
        {
            'data': {
                'id': met_pairs[0] + '_' + met_pairs[1],
                'source': met_pairs[0],
                'target': met_pairs[1],
                'label': str(np.round(met_pairs[2], 2))
            }
        })


    #*************************************Output************************************

    #*******************************************************************************

    # write to json
    print('FVA time: ' + str(end - start))
    print('Writing files...')
    directory ='data/conditions/' + dirName
    if not os.path.exists(directory):
        os.makedirs(directory)

    df = pd.DataFrame(data=reduced_met_orders)
    df.to_csv(directory + "/metaboliteOrders.csv", index=0)

    with open(directory + '/Parameters.txt', 'w') as readme:
        readme.write('Diectory: ' + dirName + '\n'
                     +'T: ' + str(T) + ' K' + '\n'
                     +'dG0 uncertainty threshold: ' + str(alpha) + '\n'
                     +'dG0 allowed error fraction: ' + str(gamma) + '\n'
                     +'Biomass threshold: ' + str(beta) + '\n'
                     +'Gibbs_eps: ' + str(Gibbs_eps) + ' kJ/mol' + '\n'
                     +'X_min: ' + str(x_min) + ' M' + '\n'
                     +'X_max: ' + str(x_max) + ' M' + '\n'
                     +'Notes: ' + notes
                    )
    with open('data/conditions/' + dirName + '/' + fileName, 'w') as outfile:
        outfile.write("data = ")
        json.dump(data, outfile)

    with open('data/graphData.json', 'w') as graphdata:
        graphdata.write("data = ")
        json.dump(data, graphdata)
