if __name__ == '__main__':

    # import numpy as np
    # import pandas as pd
    # import json
    # import time
    # import os
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
    # import networkx as nx

    # from parameters import *
    # import data
    #from classifyReactions import *
    #from MIPmodel import *
    #from candidatePairs import *


    # GEM = data.prepareGEM(work_directory=work_directory,
    #                 carbon_source=carbon_source,
    #                 uptake_rate=uptake_rate, biomass_threshold=beta,
    #                 loopless_fva=loopless_fva,
    #                 biomass_reaction_id='BIOMASS_Ec_iJO1366_core_53p95M')
    #
    # dG0_data = data.getFreeEnergyData(GEM,
    #                      work_directory=work_directory,
    #                      pH_i=pH_i,
    #                      Ionic_strength=Ionic_strength,
    #                      dG0_uncertainty_threshold=alpha)

    #from candidatePairs import *
    #import metaboliteOrders
    import output

    # Got stuck on pair 560 / 4400. I should stop a MILP if time exceeds 1 min or so
    # and discard this pair from the analysis. Also, why not write a wrapper function
    # to transform numpy matrix formulation into gurobipy input? I could use it in the future!
#
# import numpy as np
# from MIPmodel import *
# from gurobipyModelGenerator import constructGurobiModel, updateObjective
#
# c = np.zeros(N_mets + N_rxns_dG0 + N_unfixed_rxns_dG0 + N_rxns)
# c[[2, 23]] = [1, -1]
# G = np.array(G)
# h = np.array(h).flatten()
# A = np.array(A)
# b = np.array(b).flatten()
# binaryVariables = np.array(BinaryVariables)
#
# model = constructGurobiModel(c, G, h, A, b, sense='min',
#                              binaryVariables=binaryVariables, variableNames=None,
#                              modelName=None, lb=None, ub=None)
#
# Const = model.getConstrs()
# R_const = Const[-N_mets:]
# dir(R_const[0])
# # model.optimize()
# # model = setObjective(model, c, sense='min')
# # model.Status # 2 is optimal
# #
# # model.objVal()
# # model.printStats()
# # model.printQuality()
# # Vars = model.getVars()
# # binaryVarValues = [var.x for var in Vars if var.VType == 'B']
# # binaryVarValues
# # print(model.x)
# R_const[0].RHS = 0
# R_const[0].setAttr('rhs', 0)
#
# R_const[0].RHS
# model.update()
#
# x = model.getConstrs()[-N_mets:]
# x[0].RHS
