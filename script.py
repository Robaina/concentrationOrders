if __name__ == '__main__':

    from equilibrator_api import ComponentContribution, Reaction, ReactionMatcher
    import numpy as np
    import pandas as pd
    import cobra
    from cobra.flux_analysis.variability import flux_variability_analysis
    import cvxopt
    import networkx as nx
    import json
    import time
    import os


    # Constants (dG0 in kJ/mol)
    R = 8.3144598 * 1e-3 # kJ. K^-1. mol^-1
    T = 310.15 # 298.15 # K

    # Parameters
    uncertainty_threshold = alpha = 2
    biomass_threshold = beta = 1
    dG0_error_fraction = gamma = 1
    Gibbs_eps = 1e-6 # kJ/mol
    x_min, x_max = 1e-7, 1.5e-2 # M (in Teppe et al 2013* they use 1e-5, 1e-1!!
	#https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0075370)
    fileName = 'graphData.json'
    dirName = 'Cond5'
    notes = "Original iJO1366 flux bounds, added maximum ratio, new formulation"


    #******************************* Load model, dG0********************************

    #*******************************************************************************

    # Load equilibrator data
    eq_api = ComponentContribution(pH=7.5, ionic_strength=0.25) 

    # Get iJO1366 reactions as KEGG reaction strings
    iJO1366ToKEGG = pd.read_csv(
        'equilibrator-api/src/equilibrator_api/data/iJO1366_reactions.csv',
                               index_col=0)

    iJO1366 = cobra.io.load_json_model('iJO1366.json')
    # Make internal reactions reversible

    # Remove blocked reactions
    blockedRxns = cobra.flux_analysis.find_blocked_reactions(
        iJO1366, zero_cutoff=1e-9, open_exchanges=False)

    for rxn in blockedRxns:
        iJO1366.reactions.get_by_id(rxn).remove_from_model(remove_orphans=True)

    iJO1366_rxns = [rxn.id.lower() for rxn in iJO1366.reactions]


    # Prepare dictionary with dG0
    rxnGibbs = {}
    for rxn_id in iJO1366_rxns:
        try:
            rxn = Reaction.parse_formula(iJO1366ToKEGG.loc[rxn_id].item())
            dG0_prime, dG0_uncertainty = eq_api.dG0_prime(rxn)
            rxnGibbs[rxn_id] = [dG0_prime, dG0_uncertainty]
        except Exception:
            pass

    # Run fva
    print('Running flux variability analysis...')
    start = time.time()
    fva = flux_variability_analysis(iJO1366, reaction_list=None,
                                  fraction_of_optimum=beta, processes=2,
                                  loopless=True).round(decimals=4)

    end = time.time()
    fva.to_csv('iJO1366fva.csv')

    # Load pre-computed fva solution
    # fva = pd.read_csv('iJO1366fva.csv', index_col=0)


    #**************************** Irreversible reactions****************************
    #*******************************************************************************  

    # loop over fva limits
    Irr_rxns = {}
    for rxn_id in fva.index:
        v_min, v_max = fva.loc[rxn_id].minimum, fva.loc[rxn_id].maximum
        if v_min < 0 and v_max <= 0:
            Irr_rxns[rxn_id.lower()] = 'backward'
        elif v_min >= 0 and v_max > 0:
            Irr_rxns[rxn_id.lower()] = 'forward'


    # Get stoichiometric matrix (irreversible reactions with Gibb's energy data)
    GEM = iJO1366.copy()

    def notIrreversibleWithGibbsEnergy(rxn):
        rxn_id = rxn.id.lower()
        cond = ( (rxn_id not in Irr_rxns.keys())
                or (rxn_id in Irr_rxns.keys() and rxn_id not in rxnGibbs.keys())
                or (rxn_id in rxnGibbs.keys() 
                    and rxnGibbs[rxn_id][1] > alpha * rxnGibbs[rxn_id][0]) )
        return cond

    for rxn in iJO1366.reactions:
        if notIrreversibleWithGibbsEnergy(rxn):
            GEM.reactions.get_by_id(rxn.id).remove_from_model(remove_orphans=False)

    N = cobra.util.array.create_stoichiometric_matrix(GEM, array_type='dense')
    Irr_rxns_with_dG0 = [rxn.id.lower() for rxn in GEM.reactions]
    N_Irr = pd.DataFrame(data=N, index=[met.id for met in GEM.metabolites],
                         columns=Irr_rxns_with_dG0)

    # Remove orphan metabolites (no reaction)
    for met in N_Irr.index:
        if sum(abs(N_Irr.loc[met])) == 0:
            N_Irr = N_Irr.drop(met, axis=0)


    # Change direction of backward reactions
    for rxn in N_Irr.columns:
        if Irr_rxns[rxn] == 'backward':
            N_Irr[rxn] *= -1
            N_Irr.rename(columns={rxn: rxn + '_back'}, inplace=True)

    N_met_pairs = int(0.5 * len(N_Irr.index) * (len(N_Irr.index) - 1))
    print('There are ' + str(N_met_pairs) + ' metabolite pairs')


    #**************************** Linear program************************************

    #*******************************************************************************

    # Standard in cvxopt is Ax <= b so have to change signs and add epsilon to rhs
    n_rxns, n_mets = N_Irr.values.transpose().shape

    # Bounds
    logx_min = (np.log(x_min)) * np.ones((n_mets, 1))
    logx_max = (np.log(x_max)) * np.ones((n_mets, 1))

    # Construct vectors: dG0min, dG0max
    dG0min, dG0max = np.zeros((n_rxns, 1)), np.zeros((n_rxns, 1))
    for i, rxn in enumerate(Irr_rxns_with_dG0):
        dG0_i, dG0_error_i = rxnGibbs[rxn]
        dG0min[i] = dG0_i - gamma * dG0_error_i
        dG0max[i] = dG0_i + gamma * dG0_error_i

    A0 = np.hstack((N_Irr.values.transpose(), (1 / (R*T)) * np.identity(n_rxns)))
    A1 = np.hstack((np.zeros((n_rxns, n_mets)), -np.identity(n_rxns)))
    A2 = np.hstack((np.zeros((n_rxns, n_mets)), np.identity(n_rxns)))
    A3 = np.hstack((-np.identity(n_mets), np.zeros((n_mets, n_rxns))))
    A4 = np.hstack((np.identity(n_mets), np.zeros((n_mets, n_rxns))))

    A = cvxopt.matrix(np.vstack((A0, A1, A2, A3, A4)))
    b = cvxopt.matrix(np.vstack((-Gibbs_eps * np.ones((n_rxns, 1)),
                                 -dG0min, dG0max,
                                 -logx_min, logx_max)))

    # Start iteration over metabolite pairs (variable components)
    met_orders = []
    for p in range(n_mets):
        for q in range(n_mets):
            if p != q:
                # Find minimum ratio
                c = np.zeros(n_mets + n_rxns)
                c[[p, q]] = [1, -1]
                c = cvxopt.matrix(c)

                res = cvxopt.solvers.lp(c, A, b, solver='glpk',
                                        options={'glpk':{'msg_lev':'GLP_MSG_OFF'}})
                z = res['primal objective']
                if z > 0:
                    met_i, met_j = N_Irr.index[[p, q]]
                    min_ratio = np.e**z

                    # Find maximum ratio
                    c *= -1
                    res = cvxopt.solvers.lp(c, A, b, solver='glpk',
                                            options={'glpk':{'msg_lev':'GLP_MSG_OFF'}})

                    z = -res['primal objective']
                    max_ratio = np.e**z

                    met_orders.append([met_i, met_j, min_ratio, max_ratio])

					
    #************************************* Graph************************************

    #*******************************************************************************

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