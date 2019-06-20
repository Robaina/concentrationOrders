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


from parameters import *
import data

GEM = data.prepareGEM(work_directory=work_directory,
                carbon_source=carbon_source,
                uptake_rate=uptake_rate, biomass_threshold=beta,
                loopless_fva=loopless_fva,
                biomass_reaction_id='BIOMASS_Ec_iJO1366_core_53p95M')

GEM_rxns = np.array([rxn.id for rxn in GEM.reactions])
GEM_mets = np.array([met.id for met in GEM.metabolites])

dG0_data = data.getFreeEnergyData(GEM,
                     work_directory=work_directory,
                     pH_i=pH_i,
                     Ionic_strength=Ionic_strength,
                     dG0_uncertainty_threshold=alpha)

# Dependencies: GEM

########################## Categorize reactions###############################
def convertRxnIDtoIndex(GEM, rxn_id):
    return [rxn.id for rxn in GEM.reactions].index(rxn_id)
def convertMetIDtoIndex(GEM, met_id):
    return [met.id for met in GEM.metabolites].index(met_id)

Fixed_rxns_with_dG0 = [rxn.id for rxn in GEM.reactions
                  if (rxn.lower_bound > 0 and rxn.id in dG0_data.keys())]
Irr_rxns_with_dG0 = [rxn.id for rxn in GEM.reactions
                     if ('forward' not in rxn.id and 'backward' not in rxn.id
                         and rxn.id in dG0_data.keys()
                         and rxn.id not in Fixed_rxns_with_dG0)]
For_rxns_with_dG0 = [rxn.id for rxn in GEM.reactions
                     if ('forward' in rxn.id and rxn.id in dG0_data.keys()
                         and rxn.id not in Fixed_rxns_with_dG0)]
Back_rxns_with_dG0 = [rxn.id for rxn in GEM.reactions
                     if ('backward' in rxn.id and rxn.id in dG0_data.keys()
                         and rxn.id not in Fixed_rxns_with_dG0)]

Unfixed_rxns_with_dG0 = Irr_rxns_with_dG0 + For_rxns_with_dG0 + Back_rxns_with_dG0

Rest_rxns = [rxn.id for rxn in GEM.reactions
             if (rxn.id not in Unfixed_rxns_with_dG0
                 and rxn.id not in Fixed_rxns_with_dG0)]

Fixed_rxns_with_dG0_Idx = [convertRxnIDtoIndex(GEM, id) for id in Fixed_rxns_with_dG0]
Unfixed_rxns_with_dG0_Idx = [convertRxnIDtoIndex(GEM, id) for id in Unfixed_rxns_with_dG0]
Irr_rxns_with_dG0_Idx = [convertRxnIDtoIndex(GEM, id) for id in Irr_rxns_with_dG0]
For_rxns_with_dG0_Idx = [convertRxnIDtoIndex(GEM, id) for id in For_rxns_with_dG0]
Back_rxns_with_dG0_Idx = [convertRxnIDtoIndex(GEM, id) for id in Back_rxns_with_dG0]
Rest_rxns_Idx = [convertRxnIDtoIndex(GEM, id) for id in Rest_rxns]

# Extract stoichiometric matrix
N_original_order = cobra.util.array.create_stoichiometric_matrix(GEM, array_type='dense')
N_mets, N_rxns = N_original_order.shape

# Re-order reactions: Fixed_rxns_with_dG0, Irr_rxns_with_dG0, For_rxns_with_dG0,
# Back_rxns_with_dG0, Rest_rxns
re_ordered_rxns = np.concatenate((Fixed_rxns_with_dG0, Irr_rxns_with_dG0, For_rxns_with_dG0,
                                     Back_rxns_with_dG0, Rest_rxns))
re_ordered_indices = [convertRxnIDtoIndex(GEM, id) for id in re_ordered_rxns]
rxns_with_dG0 = np.concatenate((Fixed_rxns_with_dG0, Irr_rxns_with_dG0, For_rxns_with_dG0,
                                     Back_rxns_with_dG0))
rxns_with_dG0_indices = [convertRxnIDtoIndex(GEM, id) for id in rxns_with_dG0]

N = N_original_order[:, re_ordered_indices]

# Metabolites which participate in reactions with available dG0
dG0_constrained_mets = [met.id for met in GEM.metabolites
                        if any(rxn_id in rxns_with_dG0.tolist()
                               for rxn_id in [rxn.id for rxn in met.reactions])]
N_met_pairs = int(0.5 * len(dG0_constrained_mets) * (len(dG0_constrained_mets) - 1))

# # Define dimensions of submatrices
# N_rev_with_dG0 = len(For_rxns_with_dG0)
# N_rxns_no_dG0 = len(Rest_rxns)
# N_irr_with_dG0 = len(Irr_rxns_with_dG0)
# N_fixed_rxns_dG0 = len(Fixed_rxns_with_dG0)
# N_unfixed_rxns_dG0 = len(Unfixed_rxns_with_dG0)
#
# print(('There are ' + str(N_rxns - N_rxns_no_dG0)
#       + ' reactions with dG0 data and '
#       + str(N_met_pairs) + ' metabolite pairs and ' + str(N_fixed_rxns_dG0)
#       + ' fixed reactions with dG0'))

# Write file with fixed reactions with dG0
# data = [id.lower() for id in Fixed_rxns_with_dG0]
# np.savetxt('Fixed_Irr_dG0_MILP.csv', data, delimiter=',', fmt='%s')
