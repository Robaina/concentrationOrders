import numpy as np
import pandas as pd
import json
import time
import os
import networkx as nx
import parameters as par


def buildDAGadjacencyMatrix(GEM, met_orders):
    # Build full graph adjacency matrix
    total_mets = np.unique(
        np.array(
            [[met_pairs[0], met_pairs[1]]
             for met_pairs in met_orders]).flatten()).tolist()
    N_nodes = len(total_mets)

    A = np.zeros((N_nodes, N_nodes))
    for met_pair in met_orders:
        i, j = total_mets.index(
            met_pair[0]), total_mets.index(met_pair[1])
        A[i, j] = met_pair[2]

    # Get transitive reduction of the graph
    G = nx.from_numpy_matrix(A, create_using=nx.DiGraph())
    G_plus = nx.transitive_reduction(G)

    # Get adjacency matrix of the transitive reduction
    A_plus = nx.to_numpy_matrix(G)
    return (A_plus, total_mets)


def removeEdgesToProton(A_plus, total_mets):
    """Remove edges between h_c and a metabolite
    if this is the only edge of this metabolite to
    avoid saturating the graph.
    """
    def onlyProtonMetaboliteEdge(met_i, met_j):
        onlyToProton = ((total_mets[met_i] == 'h_c' and sum(A_plus[:, met_j]) == 0)
                        or (total_mets[met_j] == 'h_c' and sum(A_plus[:, met_i]) == 0))
        return onlyToProton
    
    def edgeToProton(met_i, met_j):
        return ((total_mets[met_i] == 'h_c' and total_mets[met_j] != 'h_e')
                or (total_mets[met_j] == 'h_e' and total_mets[met_j] != 'h_c'))
    
    N_nodes = A_plus.shape[0]
    reduced_met_orders = []
    for i in range(N_nodes):
        for j in range(N_nodes):
#             if not onlyProtonMetaboliteEdge(i, j):
#             if not edgeToProton(i, j):
            if total_mets[i] != 'h_c' and total_mets[j] != 'h_c':
               min_ratio = A_plus[i, j]
#             if min_ratio > 1:
               if A_plus[i, j] != 0:
                   reduced_met_orders.append(
                                [total_mets[i], total_mets[j], min_ratio])
    return reduced_met_orders


def buildGraphJSONData(GEM, ordered_met_pairs):
    # Build data dictionary (using reduced_met_orders)
    data = {}
    data['nodes'] = []
    data['edges'] = []
    
    total_mets = np.unique(
        np.array(
            [[met_pairs[0], met_pairs[1]]
             for met_pairs in ordered_met_pairs]).flatten()).tolist()

    for met in total_mets:
        data['nodes'].append(
        {
            'data': {
                'id': met,
                'label': met,
                'name': GEM.metabolites.get_by_id(met).name
            }
        })

    for met_pairs in ordered_met_pairs:
        data['edges'].append(
        {
            'data': {
                'id': met_pairs[0] + '_' + met_pairs[1],
                'source': met_pairs[0],
                'target': met_pairs[1],
                'label': str(np.round(met_pairs[2], 3))
            }
        })
    return data


def writeParametersLog(data):
     with open(f'{par.work_directory}/{par.directory}/Parameters.txt', 'w') as file:
        file.write(time.strftime("%a, %d %b %Y %H:%M:%S +0000", time.gmtime()) + '\n'
                     + 'Condition name: ' + par.dirName + '\n'
                     + 'Number of ordered pairs: {}'.format(len(data['edges'])) + '\n'
                     + ' ' + '\n'
                     + 'Parameter values:' + '\n'
                     + '-----------------' + '\n'
                     + 'T: ' + str(par.T) + ' K' + '\n'
                     + 'Intracellular pH: ' + str(par.pH_i) + ' ± ' + str(par.pH_deviation) + '\n'
                     + 'Maximum internal O_2: ' + str(par.maxo2) + ' mM \n'
                     + 'dG0 uncertainty threshold: ' + str(par.alpha) + '\n'
                     + 'Carbon uptake reaction: ' + par.carbon_source + '\n'
                     + 'Uptake rate: ' + str(par.uptake_rate) + ' mmol.min^-1.gDW^-1' + '\n'
                     + 'Biomass threshold: ' + str(par.beta) + '\n'
                     + 'dG0_eps: ' + str(-par.dG0_eps) + ' kJ/mmol' + '\n'
                     + 'X_min: ' + str(par.x_min) + ' mM' + '\n'
                     + 'X_max: ' + str(par.x_max) + ' mM' + '\n'
                     + ' ' + '\n'
                     + 'Notes: ' + '\n'
                     + '------' + '\n'
                     + par.notes)


def writeOutputFiles(data):
    print(' ')
    print('Writing files...')
    if not os.path.exists(par.work_directory + '/' + par.directory):
        os.makedirs(par.work_directory + '/' + par.directory)

    writeParametersLog(data)

    with open((par.work_directory + '/data/conditions/' 
               + par.dirName + '/' + par.fileName), 'w') as outfile:
        outfile.write("data = ")
        json.dump(data, outfile)

    with open(par.work_directory + '/data/' + par.fileName, 'w') as graphdata:
        graphdata.write("data = ")
        json.dump(data, graphdata)
