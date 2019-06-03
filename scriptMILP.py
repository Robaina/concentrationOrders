if __name__ == '__main__':

    from equilibrator_api import ComponentContribution, Reaction
    import numpy as np
    import pandas as pd
    import cobra
    from cobra.core.reaction import Reaction as cobraReaction
    from cobra.flux_analysis.variability import flux_variability_analysis
    from cobra.util.solver import set_objective
    import cvxopt
    from cvxopt import glpk
    import networkx as nx
    import json
    import time
    import os
    from six import iteritems
    import re

    # Constants (dG0 in kJ/mol) original dG0 data units
    R = 8.3144598 * 1e-6 # kJ. K^-1. mmol^-1
    T = 310.15 # 298.15 # K

    # Parameters
    loopless_fva = False
    numerator_metabolites = p_mets = None
    denominator_metabolites = q_mets = None
    dG0_uncertainty_threshold = alpha = 0.1 # filter noisy dG0 data, percentace of mean
    carbon_source = 'EX_glc__D_e'
    uptake_rate = 10 # mmol.min^-1.gDW^-1
    biomass_threshold = beta = 0.95
    dG0_error_fraction = gamma = 1 # allowed dG0 deviation from measured values
    Gibbs_eps = 1e-8 # kJ/mmol
    M = 1e8 # big number to unconstrain reactions with no dG0 data
    x_min, x_max = 1e-4, 2e2 # mM
    Intracellular_pH = pH_i = 7.3
    pH_deviation = delta_pH = 0.3
    maximum_internal_O2 = maxo2 = 2e2#9e-3 # mM
    # Fluxes in mmol.gDW^-1.min^-1
    # (in Teppe et al 2013* they use 1e-5, 1e-1!
    # https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0075370)
    fileName = 'graphData.json'
    dirName = 'Cond1'
    directory ='data/conditions/' + dirName
    notes = 'Growth on Glucose, MILP formulation!'

    # ******************************* Load model, dG0********************************
    # *******************************************************************************

    def updateExchangeReactionBounds(GEM, uptakeRxn='EX_glc__D_e', carbonUptakeRate=10):
        """
        Update exchange reaction bounds to simulate appropriate growth medium
        conditions.
        """
        def isOrganicExchange(ID):
            compoundAtoms = list(GEM.reactions.get_by_id(ID).reactants[0].formula)
            cond = (('C' in compoundAtoms)
                    & ('H' in compoundAtoms)
                    & ('o' not in compoundAtoms))  # discards cobalamine
            return cond

        ExchangeRxnIDs = [rxn.id for rxn in GEM.exchanges if len(rxn.products) == 0]
        for ID in ExchangeRxnIDs:
            try:
                if isOrganicExchange(ID):
                    GEM.reactions.get_by_id(ID).lower_bound = 0
            except Exception:
                pass
        GEM.reactions.get_by_id(uptakeRxn).lower_bound = -carbonUptakeRate
        # Temporarily remove constraint on ATPM
        GEM.reactions.get_by_id('ATPM').lower_bound = 0

    def convert_to_irreversible_SRE(cobra_model):
        """Split reversible reactions into two irreversible reactions: one going in
        the forward direction, the other in the backward direction. In this manner,
        all reactions in the model carry non-negative flux values. Forward reactions
        are tagged as "forward" while backward reactions as "backward".

        cobra_model: A Model object which will be modified in place.
        Modified from the deprecated cobrapy version by Semidan Robaina,
        February 2019.
        """
        reactions_to_add = []
        coefficients = {}
        def onlyBackward(reaction):
            return reaction.lower_bound < 0 and reaction.upper_bound <= 0

        def backwardAndForward(reaction):
            return reaction.lower_bound < 0 and reaction.upper_bound > 0

        def changeDirection(reaction):
            def swapSign(number):
                return -number
            lb = swapSign(reaction.upper_bound)
            ub = swapSign(reaction.lower_bound)
            reaction.lower_bound = lb
            reaction.upper_bound = ub
            reaction.objective_coefficient * -1
            reaction.notes["reflection"] = 'reversed, only backward'
            for met in reaction._metabolites.keys():
                reaction._metabolites[met] *= -1

        def createBackwardReaction(reaction):
            backward_reaction = cobraReaction(reaction.id + '_backward')
            backward_reaction.lower_bound = 0
            backward_reaction.upper_bound = -reaction.lower_bound
            reaction_dict = {k: v * -1
                             for k, v in iteritems(reaction._metabolites)}
            backward_reaction.add_metabolites(reaction_dict)
            backward_reaction._model = reaction._model
            backward_reaction._genes = reaction._genes
            for gene in reaction._genes:
                gene._reaction.add(backward_reaction)
            backward_reaction.subsystem = reaction.subsystem
            backward_reaction.name = reaction.name + '_backward'
            backward_reaction._gene_reaction_rule = reaction._gene_reaction_rule
            coefficients[backward_reaction] = reaction.objective_coefficient * -1
            return backward_reaction

        for reaction in cobra_model.reactions:
            if onlyBackward(reaction):
                changeDirection(reaction)
            elif backwardAndForward(reaction):
                backward_reaction = createBackwardReaction(reaction)
                reactions_to_add.append(backward_reaction)
                reaction.id += '_forward'
                reaction.name += '_forward'
                reaction.lower_bound = 0

        cobra_model.add_reactions(reactions_to_add)
        set_objective(cobra_model, coefficients, additive=True)

    # Load model (Update exchange reaction bounds)
    GEM = cobra.io.load_json_model('iJO1366.json')
    updateExchangeReactionBounds(GEM, uptakeRxn=carbon_source, carbonUptakeRate=uptake_rate)

    # Find non-default irreversible reactions
    print('Running flux variability analysis...')
    fva = flux_variability_analysis(GEM, reaction_list=None,
                                    fraction_of_optimum=beta, processes=2,
                                    loopless=loopless_fva).round(decimals=8)
    # fva.to_csv('iJO1366fva.csv')
    GEM.optimize()
    biomass_reaction = GEM.reactions.get_by_id('BIOMASS_Ec_iJO1366_core_53p95M')
    biomass_reaction.lower_bound = beta * GEM.objective.value

    # Update lower bound of fixed reactions and eliminate blocked ones
    for rxn_id in fva.index:
        v_min, v_max = fva.loc[rxn_id].minimum, fva.loc[rxn_id].maximum
        if v_min > 0 or v_max < 0:
            GEM.reactions.get_by_id(rxn_id).lower_bound = v_min
            GEM.reactions.get_by_id(rxn_id).upper_bound = v_max
        if v_min == 0 and v_max == 0:
            GEM.reactions.get_by_id(rxn_id).remove_from_model(remove_orphans=True)

    # Convert model to irreversible
    convert_to_irreversible_SRE(GEM)

    # Remove unused metabolites
    _ = cobra.manipulation.delete.prune_unused_metabolites(GEM)

    # Load equilibrator data
    eq_api = ComponentContribution(pH=pH_i, ionic_strength=0.25)

    # Get iJO1366 reactions as KEGG reaction strings
    GEM2KEGG = pd.read_csv(
        'equilibrator-api/src/equilibrator_api/data/iJO1366_reactions.csv',
                               index_col=0)

    # Prepare dictionary with dG0
    def originallyIrreversible(reaction_id):
        return 'forward' not in reaction_id and 'backward' not in reaction_id

    dG0_data = {}
    GEM_rxns = np.array([rxn.id for rxn in GEM.reactions])
    GEM_mets = np.array([met.id for met in GEM.metabolites])

    for rxn_id in GEM_rxns:
        try:
            if originallyIrreversible(rxn_id):
                id, direction = rxn_id, ''
            else:
                id, direction = re.split('_(forward|backward)', rxn_id)[:2]

            rxn = Reaction.parse_formula(GEM2KEGG.loc[id.lower()].item())
            dG0_prime, dG0_uncertainty = np.array(eq_api.dG0_prime(rxn)) * 1e-3 # convert to kJ/mmol
            if dG0_uncertainty < abs(alpha * dG0_prime): # remove uncertain dG0 data
                if 'backward' in direction:
                    dG0_data[rxn_id] = [-dG0_prime, dG0_uncertainty]
                else:
                    dG0_data[rxn_id] = [dG0_prime, dG0_uncertainty]
        except Exception:
            pass

    # ************************************Prepare MILP************************************
    # ************************************************************************************

    # Categorize reactions
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

    # Write file with fixed reactions with dG0
    # data = [id.lower() for id in Fixed_rxns_with_dG0]
    # np.savetxt('Fixed_Irr_dG0_MILP.csv', data, delimiter=',', fmt='%s')

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
    BinaryVariables = set(range(N_mets + N_rxns_dG0, N_mets + N_rxns_dG0 + N_unfixed_rxns_dG0))

    # ******************************Loop over metabolites*********************************
    # ************************************************************************************

    # Start iteration over metabolite pairs (variable components)
    met_orders = []
    if p_mets is None:
        p_mets = [convertMetIDtoIndex(GEM, met_id) for met_id in dG0_constrained_mets]
    else:
        p_mets = [convertMetIDtoIndex(GEM, met_id) for met_id in p_mets]
    if q_mets is None:
        q_mets = [convertMetIDtoIndex(GEM, met_id) for met_id in dG0_constrained_mets]
    else:
        q_mets = [convertMetIDtoIndex(GEM, met_id) for met_id in q_mets]

    start = time.time()
    print('Iterating over metabolite pairs...')
    for p in p_mets:
        for q in q_mets:
            if p != q:
                # Find minimum ratio
                c = np.zeros(N_mets + N_rxns_dG0 + N_unfixed_rxns_dG0 + N_rxns)
                c[[p, q]] = [1, -1]
                c = cvxopt.matrix(c)

                res = glpk.ilp(c, G, h, A, b, B=set(BinaryVariables),
                               options={'msg_lev':'GLP_MSG_OFF'})
                x = res[1]
                z = x[p] - x[q]
                if z > 0:
                    met_i, met_j = GEM_mets[[p, q]]
                    min_ratio = np.e**z

                    # Find maximum ratio
                    c *= -1
                    res = glpk.ilp(c, G, h, A, b, B=set(BinaryVariables),
                                   options={'msg_lev':'GLP_MSG_OFF'})

                    x = res[1]
                    z = x[p] - x[q]
                    max_ratio = np.e**z
                    # added 27 Frebruary
                    min_ratio_j_i = np.e**(x[q] - x[p])

                    met_orders.append([met_i, met_j, min_ratio, max_ratio, min_ratio_j_i])

    end = time.time()
    df = pd.DataFrame(data=met_orders)
    df.to_csv(directory + "/FullmetaboliteOrders.csv", index=0)
    # ************************************* Graph*****************************************
    # ************************************************************************************

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
    print('Total time: ' + str(end - start) + ' seconds')
    print('Writing files...')
    
    if not os.path.exists(directory):
        os.makedirs(directory)

    df = pd.DataFrame(data=reduced_met_orders)
    df.to_csv(directory + "/metaboliteOrders.csv", index=0)

    with open(directory + '/Parameters.txt', 'w') as readme:
        readme.write('Directory: ' + dirName + '\n'
                     +time.strftime("%a, %d %b %Y %H:%M:%S +0000", time.gmtime()) + '\n'
                     +'T: ' + str(T) + ' K' + '\n'
					 +'Intracellular pH: ' + str(pH_i) + ' ± ' + str(pH_deviation) + '\n'
                     +'Maximum internal O_2: ' + str(maxo2) + ' mM \n'
                     +'dG0 uncertainty threshold: ' + str(alpha) + '\n'
                     +'dG0 allowed error fraction: ' + str(gamma) + '\n'
                     +'Carbon uptake reaction: ' + carbon_source + '\n'
                     +'Uptake rate: ' + str(uptake_rate) + ' mmol.min^-1.gDW^-1' + '\n'
                     +'Biomass threshold: ' + str(beta) + '\n'
                     +'Gibbs_eps: ' + str(Gibbs_eps) + ' kJ/mmol' + '\n'
                     +'X_min: ' + str(x_min) + ' mM' + '\n'
                     +'X_max: ' + str(x_max) + ' mM' + '\n'
                     +'Loopless FVA: ' +  str(loopless_fva) + '\n'
                     +'Notes: ' + notes
                    )
    with open('data/conditions/' + dirName + '/' + fileName, 'w') as outfile:
        outfile.write("data = ")
        json.dump(data, outfile)

    with open('data/graphData.json', 'w') as graphdata:
        graphdata.write("data = ")
        json.dump(data, graphdata)


# Testing:********************************************************************************
# ****************************************************************************************
# res
# print(met_orders)
#
# def getMILPIndexFromReactionID(rxn_id):
#     return N_mets + N_rxns_dG0 + N_unfixed_rxns_dG0 + re_ordered_rxns.tolist().index(rxn_id)
# v = [x[getMILPIndexFromReactionID(id)] for id in GEM_rxns]
# sum(np.dot(N_original_order, v))
# 'ACCOAC', 'ATPM', 'BIOMASS_Ec_iJO1366_core_53p95M'
# x[getMILPIndexFromReactionID('ATPM')]
# import numpy as np
# MILP_data = np.genfromtxt('Fixed_Irr_dG0_MILP.csv', delimiter=',', dtype='U20')
# LP_data = np.genfromtxt('Fixed_Irr_dG0_LP.csv', delimiter=',', dtype='U20')
# unique_LP = np.setdiff1d(LP_data, MILP_data)
# GEM.reactions.get_by_id('FBA')
#
# 'ACCOAL' in Irr_rxns_with_dG0
