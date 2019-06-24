import numpy as np
import pandas as pd
# import json
import time
# import os
# from six import iteritems
# import re

# import cobra
# from cobra.core.reaction import Reaction as cobraReaction
# from cobra.flux_analysis.variability import flux_variability_analysis
# from cobra.util.solver import set_objective
# from equilibrator_api import ComponentContribution, Reaction
import cvxopt
from cvxopt import glpk
# import networkx as nx

from MIPmodel import *
from candidatePairs import candidate_pairs
from gurobipyModelGenerator import constructGurobiModel, updateObjective

met_orders = []
start = time.time()
print(' ')
print('Finding metabolite pairs...')
n = 0
N_pairs = len(candidate_pairs)


# Construct gurobi model
c = np.zeros(N_mets + N_rxns_dG0 + N_unfixed_rxns_dG0 + N_rxns)
c[[1, 2]] = [1, -1]
G = np.array(G)
h = np.array(h).flatten()
A = np.array(A)
b = np.array(b).flatten()
binaryVariables = np.array(BinaryVariables)

model = constructGurobiModel(c, G, h, A, b, sense='min',
                             binaryVariables=binaryVariables, variableNames=None,
                             modelName=None, lb=None, ub=None)
model.setParam('OutputFlag', False)


for pair in candidate_pairs:
    p, q = pair

    # Find minimum ratio
    c = np.zeros(N_mets + N_rxns_dG0 + N_unfixed_rxns_dG0 + N_rxns)
    c[[p, q]] = [1, -1]
    # c = cvxopt.matrix(c)

    # res = glpk.ilp(c, G, h, A, b, B=set(BinaryVariables),
    #                options={'msg_lev': 'GLP_MSG_OFF', 'tm_lim': 5000})
    # (status, x) = res

    model = updateObjective(model, c)
    model.optimize()
    x, status = model.X, model.Status

    z = x[p] - x[q]

    if z > 0 and status == GRB.Status.OPTIMAL:  # in 'optimal':
        met_i, met_j = GEM_mets[[p, q]]
        min_ratio = np.e**z

        # # Find maximum ratio
        # c *= -1
        # res = glpk.ilp(c, G, h, A, b, B=set(BinaryVariables),
        #                options={'msg_lev':'GLP_MSG_OFF'})
        #
        # x = res[1]
        # z = x[p] - x[q]
        # max_ratio = np.e**z
        # # added 27 Frebruary
        # min_ratio_j_i = np.e**(x[q] - x[p])

        # met_orders.append([met_i, met_j, min_ratio, max_ratio, min_ratio_j_i])

        met_orders.append([met_i, met_j, min_ratio])
    print(('Iterating... ({} of ' + str(N_pairs) + ')').format(n), end='\r')
    n += 1

end = time.time()
print(' ')
print('Total time: {:.2f} seconds'.format(end - start))
print('There are a total of {} ordered pairs'.format(len(met_orders)))
