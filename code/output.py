import numpy as np
import pandas as pd
import json
import time
import os
# from six import iteritems
# import re
#
# import cobra
# from cobra.core.reaction import Reaction as cobraReaction
# from cobra.flux_analysis.variability import flux_variability_analysis
# from cobra.util.solver import set_objective
# from equilibrator_api import ComponentContribution, Reaction
# import cvxopt
# from cvxopt import glpk
import networkx as nx

from parameters import *
from metaboliteOrders import met_orders, GEM

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
    #B[i, j] = met_pair[3]

# Get transitive reduction of the graph
G = nx.from_numpy_matrix(A, create_using=nx.DiGraph())
G_plus = nx.transitive_reduction(G)

# Get adjacency matrix of the transitive reduction
A_plus = nx.to_numpy_matrix(G)


def onlyProtonMetaboliteEdge(met_i, met_j):
    """Remove edges between h_c and a metabolite
    if this is the only edge of this metabolite to
    avoid saturating the graph.
    """
    onlyToProton = ((total_mets[met_i] == 'h_c' and sum(A_plus[:, met_j]) == 0)
                    or (total_mets[met_j] == 'h_c' and sum(A_plus[:, met_i]) == 0))
    return onlyToProton


# Get list of edges
reduced_met_orders = []
for i in range(N_nodes):
    for j in range(N_nodes):
        #if not onlyProtonMetaboliteEdge(i, j):
        min_ratio = A_plus[i, j]
        if min_ratio > 1:
            #max_ratio = B[i, j]
            reduced_met_orders.append(
                            [total_mets[i], total_mets[j], min_ratio])  # , max_ratio])

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
            'name': GEM.metabolites.get_by_id(met).name
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

#*************************************Output******************************************
#*************************************************************************************

# write to json
print(' ')
print('Writing files...')
if not os.path.exists(work_directory + directory):
    os.makedirs(work_directory + directory)

df = pd.DataFrame(data=reduced_met_orders)
df.to_csv(work_directory + directory + "/metaboliteOrders.csv", index=0)

# df = pd.DataFrame(data=met_orders)
# df.to_csv(work_directory + directory
#          + "/FullmetaboliteOrders.csv", index=0)

with open(work_directory + directory + '/Parameters.txt', 'w') as readme:
    readme.write(time.strftime("%a, %d %b %Y %H:%M:%S +0000", time.gmtime()) + '\n'
                 + 'Condition name: ' + dirName + '\n'
                 + 'Number of ordered pairs: {}'.format(len(met_orders)) + '\n'
                 + ' ' + '\n'
                 + 'Parameter values:' + '\n'
                 + '-----------------' + '\n'
                 + 'T: ' + str(T) + ' K' + '\n'
                 + 'Intracellular pH: ' + str(pH_i) + ' ± ' + str(pH_deviation) + '\n'
                 + 'Maximum internal O_2: ' + str(maxo2) + ' mM \n'
                 + 'dG0 uncertainty threshold: ' + str(alpha) + '\n'
                 + 'dG0 allowed error fraction: ' + str(gamma) + '\n'
                 + 'Carbon uptake reaction: ' + carbon_source + '\n'
                 + 'Uptake rate: ' + str(uptake_rate) + ' mmol.min^-1.gDW^-1' + '\n'
                 + 'Biomass threshold: ' + str(beta) + '\n'
                 + 'Gibbs_eps: ' + str(Gibbs_eps) + ' kJ/mmol' + '\n'
                 + 'X_min: ' + str(x_min) + ' mM' + '\n'
                 + 'X_max: ' + str(x_max) + ' mM' + '\n'
                 + 'Loopless FVA: ' +  str(loopless_fva) + '\n'
                 + ' ' + '\n'
                 + 'Notes: ' + '\n'
                 + '------' + '\n'
                 + notes)

with open(work_directory + 'data/conditions/' + dirName + '/' + fileName, 'w') as outfile:
    outfile.write("data = ")
    json.dump(data, outfile)

with open(work_directory + 'data/graphData.json', 'w') as graphdata:
    graphdata.write("data = ")
    json.dump(data, graphdata)
