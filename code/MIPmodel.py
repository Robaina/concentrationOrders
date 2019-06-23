import numpy as np
import cvxopt

# from parameters import *
from classifyReactions import *

# Define dimensions of submatrices
N_rev_with_dG0 = len(For_rxns_with_dG0)
N_rxns_no_dG0 = len(Rest_rxns)
N_irr_with_dG0 = len(Irr_rxns_with_dG0)
N_fixed_rxns_dG0 = len(Fixed_rxns_with_dG0)
N_unfixed_rxns_dG0 = len(Unfixed_rxns_with_dG0)

print(('There are ' + str(N_rxns - N_rxns_no_dG0)
      + ' reactions with dG0 data and '
      + str(N_met_pairs) + ' metabolite pairs and ' + str(N_fixed_rxns_dG0)
      + ' fixed reactions with dG0'))

# Standard in cvxopt is Ax <= b so have to change signs and add epsilon to rhs
# Bounds
logx_min = (np.log(x_min)) * np.ones((N_mets, 1))
logx_max = (np.log(x_max)) * np.ones((N_mets, 1))

# Fix internal pH (change to mM)
h_index = GEM_mets.tolist().index('h_c')
logx_min[h_index] = np.log(10**(-(pH_i + delta_pH) + 3))
logx_max[h_index] = np.log(10**(-(pH_i - delta_pH) + 3))

# Upper bound to internal o2
o2_index = GEM_mets.tolist().index('o2_c')
logx_max[o2_index] = np.log(maxo2) #mM

v_min = np.array([rxn.lower_bound for rxn in GEM.reactions])
v_max = np.array([rxn.upper_bound for rxn in GEM.reactions])

# Construct vectors: dG0min, dG0max
N_rxns_dG0 = len(rxns_with_dG0)
dG0min, dG0max = np.zeros((N_rxns_dG0, 1)), np.zeros((N_rxns_dG0, 1))

for i, rxn in enumerate(rxns_with_dG0):
    dG0_i, dG0_error_i = dG0_data[rxn]
    dG0min[i] = dG0_i - gamma * dG0_error_i
    dG0max[i] = dG0_i + gamma * dG0_error_i

# Build constraints matrix (cols: N_mets + N_rxns_dG0 + N_unfixed_rxns_dG0 + N_rxns )

# Inequality constraints (logx, dG0, y_irr, y_for, y_back, v_irr_dG0, v_for_dG0, v_back_dG0, v_no_dG0)
A0 = np.hstack((R*T*N[:, :N_fixed_rxns_dG0].transpose(), np.identity(N_fixed_rxns_dG0),
                np.zeros((N_fixed_rxns_dG0, 2 * N_unfixed_rxns_dG0 + N_rxns)))) # dG0 fixed

A1 = np.hstack((R*T*N[:, N_fixed_rxns_dG0:N_rxns_dG0].transpose(),
                np.zeros((N_unfixed_rxns_dG0, N_fixed_rxns_dG0)), np.identity(N_unfixed_rxns_dG0),
                M * np.identity(N_unfixed_rxns_dG0), np.zeros((N_unfixed_rxns_dG0, N_rxns)))) # dG0 non-fixed

A2 = np.hstack((np.zeros((N_rxns_dG0, N_mets)), -np.identity(N_rxns_dG0),
                np.zeros((N_rxns_dG0, N_unfixed_rxns_dG0 + N_rxns)))) # dG0 min

A3 = np.hstack((np.zeros((N_rxns_dG0, N_mets)), np.identity(N_rxns_dG0),
                np.zeros((N_rxns_dG0, N_unfixed_rxns_dG0 + N_rxns)))) # dG0 max

A4 = np.hstack((np.zeros((N_unfixed_rxns_dG0, N_mets + N_rxns_dG0)),
                np.diag(v_min[Unfixed_rxns_with_dG0_Idx]),
                np.zeros((N_unfixed_rxns_dG0, N_fixed_rxns_dG0)), # v min with data
                -np.identity(N_unfixed_rxns_dG0), np.zeros((N_unfixed_rxns_dG0, N_rxns_no_dG0))))

A5 = np.hstack((np.zeros((N_unfixed_rxns_dG0, N_mets + N_rxns_dG0)),
                np.diag(-v_max[Unfixed_rxns_with_dG0_Idx]),
                np.zeros((N_unfixed_rxns_dG0, N_fixed_rxns_dG0)), # v max with data
                np.identity(N_unfixed_rxns_dG0), np.zeros((N_unfixed_rxns_dG0, N_rxns_no_dG0))))

A6 = np.hstack((np.zeros((N_rev_with_dG0, N_mets + N_rxns_dG0 + N_irr_with_dG0)),
                np.identity(N_rev_with_dG0), np.identity(N_rev_with_dG0),
                np.zeros((N_rev_with_dG0, N_rxns)))) # y+ + y- <= 1

A7 = np.hstack((-np.identity(N_mets),
                np.zeros((N_mets, N_rxns_dG0 + N_unfixed_rxns_dG0 + N_rxns)))) # minlogx

A8 = np.hstack((np.identity(N_mets),
                np.zeros((N_mets, N_rxns_dG0 + N_unfixed_rxns_dG0 + N_rxns)))) # maxlogx

A9A = np.hstack((np.zeros((N_fixed_rxns_dG0, N_mets + N_rxns_dG0 + N_unfixed_rxns_dG0)),
               -np.identity(N_fixed_rxns_dG0),
                np.zeros((N_fixed_rxns_dG0, N_unfixed_rxns_dG0 + N_rxns_no_dG0)))) #vmin fixed dG0

A9B = np.hstack((np.zeros((N_fixed_rxns_dG0, N_mets + N_rxns_dG0 + N_unfixed_rxns_dG0)),
                np.identity(N_fixed_rxns_dG0),
                np.zeros((N_fixed_rxns_dG0, N_unfixed_rxns_dG0 + N_rxns_no_dG0)))) #vmax fixed dG0

A10A = np.hstack((np.zeros((N_rxns_no_dG0, N_mets + 2*N_rxns_dG0 + N_unfixed_rxns_dG0)),
               -np.identity(N_rxns_no_dG0))) #vmin no dG0

A10B = np.hstack((np.zeros((N_rxns_no_dG0, N_mets + 2*N_rxns_dG0 + N_unfixed_rxns_dG0)),
                np.identity(N_rxns_no_dG0))) #vmax no dG0

G = cvxopt.matrix(np.vstack((A0, A1, A2, A3, A4, A5, A6, A7, A8, A9A, A9B, A10A, A10B)))
h = cvxopt.matrix(np.vstack((-Gibbs_eps * np.ones((N_fixed_rxns_dG0, 1)),
                             -Gibbs_eps * np.ones((N_unfixed_rxns_dG0, 1)) + M,
                             -dG0min, dG0max,
                             np.zeros((2 * N_unfixed_rxns_dG0, 1)),
                             np.ones((N_rev_with_dG0, 1)),
                             -logx_min, logx_max,
                             -v_min[Fixed_rxns_with_dG0_Idx].reshape((-1, 1)),
                             v_max[Fixed_rxns_with_dG0_Idx].reshape((-1, 1)),
                             -v_min[Rest_rxns_Idx].reshape((-1, 1)),
                             v_max[Rest_rxns_Idx].reshape((-1, 1)))))

# Equality constraints
A = cvxopt.matrix(np.hstack((np.zeros((N_mets, N_mets + N_rxns_dG0 + N_unfixed_rxns_dG0)), N))) # Nv = 0
b = cvxopt.matrix(np.zeros((N_mets, 1)))
BinaryVariables = list(range(N_mets + N_rxns_dG0, N_mets + N_rxns_dG0 + N_unfixed_rxns_dG0))
