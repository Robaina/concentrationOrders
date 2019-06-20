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


#from parameters import *
from MIPmodel import *

"""
First need to modify the MILP model to integrate the first norm:
Add dummy variables q+, q- of size N_mets and the constraint
q+ - q- = log_xrand - log_x, with q+, q- >= 0 (q+ - q- + log_x = log_xrand)
"""

# Inequality constraints (logx, dG0, y_irr, y_for, y_back, v_irr_dG0, v_for_dG0, v_back_dG0, v_no_dG0, q+, q-)

# Add zeros to inequality constraints matrix to match new dimensions
G_temp = np.hstack((np.array(G), np.zeros((np.size(G, 0), 2 * N_mets))))

# Add bound constraint for q+, q-
q_bounds = np.vstack((
                      np.hstack((np.zeros((N_mets, np.size(G, 1))),
                      -np.identity(N_mets), np.zeros((N_mets, N_mets)))),

                      np.hstack((np.zeros((N_mets, np.size(G, 1))),
                      np.zeros((N_mets, N_mets)), -np.identity(N_mets)))
                      ))
Gp = cvxopt.matrix(np.vstack((G_temp, q_bounds)))

# Add new terms to rhs
hp = cvxopt.matrix(np.vstack((np.array(h),
                              np.zeros((2 * N_mets, 1)))))

# Add zeros to equality constraints matrix to match new dimensions
A_temp = np.hstack((np.array(A), np.zeros((np.size(A, 0), 2 * N_mets))))

# Add constraint on q+, q_
q_constraint = np.hstack((np.identity(N_mets),
                          np.zeros((N_mets, N_rxns_dG0 + N_unfixed_rxns_dG0 + N_rxns)),
                          np.identity(N_mets),
                          -np.identity(N_mets)))
Ap = cvxopt.matrix(np.vstack((A_temp, q_constraint)))

# Create objective vector: q+ + q-
c = cvxopt.matrix(np.concatenate((
                             np.zeros(N_mets + N_rxns_dG0 + N_unfixed_rxns_dG0 + N_rxns),
                             np.ones(2 * N_mets))))


# Start sampling loop
start = time.time()
logx_samples = []
n = 0
while len(logx_samples) < number_preprocessing_samples:
    try:
        # Create random logx vector
        rand_logx = (logx_max - logx_min) * np.random.rand(N_mets, 1) + logx_min

        # Add new term to equality rhs
        bp = cvxopt.matrix(np.vstack((np.array(b), rand_logx)))
        res = glpk.ilp(c, Gp, hp, Ap, bp, B=set(BinaryVariables),
                       #options={'tm_lim': 5000, 'mipgap': 10})
                       options={'msg_lev':'GLP_MSG_OFF',
                                'mipgap': 10,
                                'tm_lim': 5000})

        logx = np.array(res[1][:N_mets]).flatten()
        logx_samples.append(logx)
        n += 1
        print(n)

    except Exception:
        print('Error occurred!')
        pass

logx_samples = np.array(logx_samples)

# Check for candidate ordered metabolite pairs
def isCandidatePair(met_i_values, met_j_values):
    return all(met_i_values >= met_j_values)

candidate_pairs = []
for met_i in range(N_mets):
    for met_j in range(N_mets):
        met_i_values, met_j_values = logx_samples[:, met_i], logx_samples[:, met_j]
        if isCandidatePair(met_i_values, met_j_values):
            candidate_pairs.append([met_i, met_j])

end = time.time()
print('Total time: ' + str(end - start) + ' seconds')
print(len(candidate_pairs))
