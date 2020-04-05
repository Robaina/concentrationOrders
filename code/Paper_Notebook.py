#!/usr/bin/env python
# coding: utf-8

if __name__ == '__main__':
    from gurobipy import GRB
    from importlib import reload
    import numpy as np
    import pandas as pd
    import os
    import pickle
    import parameters as par
    import data
    import model
    import output
    reload(par)


    def saveToPickleFile(python_object, path_to_file='object.pkl'):
        """
        Save python object to pickle file
        """
        out_file = open(path_to_file,'wb')
        pickle.dump(python_object, out_file)
        out_file.close()
        
    def readFromPickleFile(path_to_file='object.pkl'):
        """
        Load python object from pickle file.
        Returns python object.
        """
        in_file = open(path_to_file,'rb')
        python_object = pickle.load(in_file)
        return python_object
    
    def getPartition(max_step, step_size=100):
        partition = []
        steps = np.arange(0, max_step, step_size)
        for i in range(len(steps) - 1):
            partition.append([steps[i], steps[i+1]])
        partition[-1][1] = max_step
        return partition

    
    GEM = data.prepareGEM(path_to_GEM=f'{par.work_directory}/{par.model}',
                          carbon_source=par.carbon_source,
                          uptake_rate=par.uptake_rate,
                          loopless=False,
                          biomass_reaction_id=par.biomass_rxn_id)

    dG0_data = data.getFreeEnergyData(GEM,
                                      work_directory=par.work_directory,
                                      pH_i=par.pH_i,
                                      Ionic_strength=par.Ionic_strength,
                                      dG0_uncertainty_threshold=par.alpha)

    # Build gurobi TFBA model
    reload(par)
    m, variables = model.buildMILPModel(GEM, dG0_data)

    # Fix internal pH (change to mM)
    logx_h_c = variables['logx']['logx_h_c']
    logx_h_c.lb = np.log(10**(-(par.pH_i + par.delta_pH) + 3))
    logx_h_c.ub = np.log(10**(-(par.pH_i - par.delta_pH) + 3))

    # Upper bound to internal h2o (mM)
    logx_h2o_c = variables['logx']['logx_h2o_e']
    logx_h2o_c.ub = np.log(1e2) # mM

    # Upper bound to internal o2 (mM)
    logx_o2_c = variables['logx']['logx_o2_c']
    logx_o2_c.ub = np.log(par.maxo2) # mM

    # Fixed external glucose concentration
    logx_glu_e = variables['logx']['logx_glu__L_e']
    logx_glu_e.lb = np.log(20) # mM 22.2
    logx_glu_e.ub = np.log(23) # mM

    # Find biomass maximum under thermodynamic constraints and constraint biomass reaction
    m.setObjective(variables['v']['v_' + par.biomass_rxn_id], GRB.MAXIMIZE)
    m.update()
    m.setParam('OutputFlag', False)
    m.optimize()

    try:
        print(f'Maximum growth rate: {m.objval:.2f} h^-1')
        max_biomass_flux = m.getVarByName('v_' + par.biomass_rxn_id).x
        m.getVarByName('v_' + par.biomass_rxn_id).lb = par.beta * max_biomass_flux
    except:
        print('Model is infeasible!')

    # Some sanity checks of the solution
    conc_sum = 0
    for var in m.getVars():
        if 'x_' in var.varName[:2]:
            conc_sum += var.x
    total_conc_internal_mets = conc_sum - m.getVarByName('x_glu__L_e').x
    print(f'Total concentration of internal metabolites: {total_conc_internal_mets:.2f} mM')
    print(f'Total sum of metabolite concentrations: {conc_sum:.2f} mM')


    # Find candidate ordered metabolite pairs
    if not os.path.exists(par.work_directory + '/' + par.directory):
        os.makedirs(par.work_directory + '/' + par.directory)
        
    m.update()
#     X_sample, candidatePairs, sampling_model = model.findCandidatePairs(
#         m, n_samples=par.n_samples)
    
#     systems_df = pd.read_excel(f'{par.work_directory}/iML1515_subsystems.xlsx')
#     central_metabolism = ['Carbohydrate metabolism'
#                       ,'Amino acid metabolism', 'Nucleotide metabolism',
#                        'Lipid metabolism']
#     central_metabolites = data.findMetabolitesInPathways(GEM, systems_df,
#                                                      central_metabolism)
#     centralCandidatePairs = []
#     for pair in candidatePairs:
#         met_i, met_j = pair
#         if met_i in central_metabolites and met_j in central_metabolites:
#             centralCandidatePairs.append(pair)
            
    
#     saveToPickleFile(X_sample, par.work_directory + '/' + par.directory + '/X_sample.pkl')
#     saveToPickleFile(centralCandidatePairs, 
#                      par.work_directory + '/' + par.directory + '/centralCandidatePairs.pkl')
        
    _, _, sampling_model = model.findCandidatePairs(m, 2)
    centralCandidatePairs = readFromPickleFile(
            path_to_file='centralCandidatePairs.pkl')

    # Evaluate metabolite orders from candidate pairs
    data_dir = par.work_directory + '/' + par.directory + '/ordered_pairs'
    if not os.path.exists(data_dir):
        os.makedirs(data_dir)
        
    partition = getPartition(len(centralCandidatePairs), step_size=100)
    for n, range_tuple in enumerate(partition):
        print(f'Evaluating ordered pairs in range: {range_tuple}')
        
        candidates = centralCandidatePairs[range_tuple[0]:range_tuple[1]]
        ordered_pairs = model.findConcentrationOrderedPairs(sampling_model, candidates)
        saveToPickleFile(ordered_pairs, path_to_file=f'{data_dir}/ordered_pairs_{n}.pkl')




    # ## 4.4 Evaluating concentration order DAG
    # A_plus, total_mets = output.buildDAGadjacencyMatrix(GEM, centralCandidatePairs)
    # reduced_met_orders = output.removeEdgesToProton(A_plus, total_mets)
    # graphData = output.buildGraphJSONData(GEM, total_mets, reduced_met_orders)
    # output.writeOutputFiles(graphData)

