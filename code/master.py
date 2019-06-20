if __name__ == '__main__':

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

    # print(GEM.reactions[0].id)
    # print(dG0_data.keys())

    from candidatePairs import *

    #print(GEM.reactions[0].id)
    #print(dG0_data.keys())


    # TODO: Fails to import cobra before importing GEM.py, cobra is not defined in GEM.py when calling it in line 85. I guess I can only define functions in modules then? Is there another alternative?
