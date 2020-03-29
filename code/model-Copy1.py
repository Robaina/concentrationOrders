import re
import numpy as np
import cobra
from gurobipy import Model, GRB
import parameters as par


def optimizeGEM(GEM, obj_rxn):
    """
    Plain old FBA using cobra model with gurobi. I made this function because the cobra
    optimize function gives error when used on split, irreversible model
    """
    model = Model('FBA_model')
    GEM_rxn_ids = [rxn.id for rxn in GEM.reactions]
    S = cobra.util.array.create_stoichiometric_matrix(GEM, array_type='dense')

    # Add variables
    v = {}
    for rxn in GEM.reactions:
        var_id = f'v_{rxn.id}'
        x = model.addVar(lb=rxn.lower_bound, ub=rxn.upper_bound,
                         obj=0.0, vtype=GRB.CONTINUOUS, name=var_id)
        v[var_id] = x

    # Add constraints
    for i, row in enumerate(S):
        c = ''
        for j, coeff in enumerate(row):
            rxn_id = GEM_rxn_ids[j]
            if coeff != 0:
                c += f'{coeff} * v["v_{rxn_id}"] +'
        c = c[:-1]
        c += '== 0'
        model.addConstr(eval(c), f'mass_balance_{GEM.metabolites[i].id}')
        
    # Set Objective
    model.setObjective(v['v_' + obj_rxn], GRB.MAXIMIZE)
    model.optimize()
    return model
            
        
def isReversible(rxn_id):
    return 'forward' in rxn_id or 'backward' in rxn_id


def getOppositeDirection(rxn_id):
    try:
        id, direction = re.split('_(forward|backward)', rxn_id)[:2]
        if 'forward' in direction:
            return f'{id}_backward'
        else:
            return f'{id}_forward'
    except Exception:
        return rxn_id


def buildMILPModel(GEM, dG0data):

    model = Model('TFA_model')
    GEM_rxn_ids = [rxn.id for rxn in GEM.reactions]
    S = cobra.util.array.create_stoichiometric_matrix(GEM, array_type='dense')

    # Add variables
    variables = {'x': {}, 'logx': {}, 'v': {}, 'dG0': {}, 'y': {}}
    for met in GEM.metabolites:
        var_id = f'logx_{met.id}'
        x = model.addVar(lb=np.log(par.x_min), ub=np.log(par.x_max),
                         obj=0.0, vtype=GRB.CONTINUOUS, name=var_id)
        variables['logx'][var_id] = x

    for rxn in GEM.reactions:
        var_id = f'v_{rxn.id}'
        x = model.addVar(lb=rxn.lower_bound, ub=rxn.upper_bound,
                         obj=0.0, vtype=GRB.CONTINUOUS, name=var_id)
        variables['v'][var_id] = x

    for rxn_id in dG0data.keys():
        var_id = f'dG0_{rxn_id}'
        lb = dG0data[rxn_id]['dG0'] - par.gamma * dG0data[rxn_id]['error']
        ub = dG0data[rxn_id]['dG0'] + par.gamma * dG0data[rxn_id]['error']
        x = model.addVar(lb=lb, ub=ub,
                         obj=0.0, vtype=GRB.CONTINUOUS, name=var_id)
        variables['dG0'][var_id] = x

    for rxn_id in dG0data.keys():
        var_id = f'y_{rxn_id}'
        x = model.addVar(lb=0, ub=1,
                         obj=0.0, vtype=GRB.BINARY, name=var_id)
        variables['y'][var_id] = x

    # Add constraints
    for i, row in enumerate(S):
        c = ''
        for j, coeff in enumerate(row):
            rxn_id = GEM_rxn_ids[j]
            if coeff != 0:
                c += f'{coeff} * variables["v"]["v_{rxn_id}"] +'
        c = c[:-1]
        c += '== 0'
        model.addConstr(eval(c), f'mass_balance_{GEM.metabolites[i].id}')

    for rxn_id in dG0data.keys():
        sum_str = ''
        rxn = GEM.reactions.get_by_id(rxn_id)
        for met in rxn.metabolites.keys():
            coeff = rxn.metabolites[met]
            sum_str += f'{coeff} * variables["logx"]["logx_{met.id}"] +'
        sum_str = sum_str[:-1]

        c = (f'variables["dG0"]["dG0_{rxn_id}"] + par.R * par.T * ({sum_str})'
             + f' - (1 - variables["y"]["y_{rxn_id}"]) * par.M <= -par.dG0_eps')
        model.addConstr(eval(c), f'dG0_{rxn_id}')

    for rxn_id in dG0data.keys():
        c_max = f'variables["v"]["v_{rxn_id}"] - variables["y"]["y_{rxn_id}"] * 1000 <= 0'
        model.addConstr(eval(c_max), f'force_max_flux_{rxn_id}')

    for rxn_id in dG0data.keys():
        if isReversible(rxn_id):
            opposite_rxn_id = getOppositeDirection(rxn_id)
            c = f'variables["y"]["y_{rxn_id}"] + variables["y"]["y_{opposite_rxn_id}"] <= 1'
            model.addConstr(eval(c), f'y_XOR_{rxn_id}')
            
    # Add maximum concentration constraint (300 mM)
#     sum_str = ''
#     for var_id in variables['logx'].keys():
#         sum_str += f'variables["logx"]["{var_id}"] +'
#     sum_str = sum_str[:-1]
#     model.addConstr(eval(c), 'totla concentration sum')
#     c_min = f'{sum_str} >= {par.min_sum_logx}'
#     model.addConstr(eval(c_min), 'min total concentration sum')
#     c_max = f'{sum_str} <= {par.max_sum_logx}'
#     model.addConstr(eval(c_max), 'max total concentration sum')


    # Adding exponential constraint x = exp(logx), only available in Gurobi 9
    for met in GEM.metabolites:
        var_id = f'x_{met.id}'
        x = model.addVar(lb=0, ub=2 * par.x_max,
                         obj=0.0, vtype=GRB.CONTINUOUS, name=var_id)
        variables['x'][var_id] = x
    
    for met in GEM.metabolites:
        logx = variables['logx'][f'logx_{met.id}']
        x = variables['x'][f'x_{met.id}']
        model.addGenConstrExp(logx, x, name=f'x_{met.id} = exp(logx_{met.id})') 
        
    sum_str = ''
    for met in GEM.metabolites:
        sum_str += f'variables["x"]["x_{met.id}"] +'
    sum_str = sum_str[:-1]
    c_min = f'{sum_str} >= {par.min_sum_x}'
    c_max = f'{sum_str} <= {par.max_sum_x}'
    model.addConstr(eval(c_min), 'min concentration sum')
    model.addConstr(eval(c_max), 'max concentration sum')
    
    
    model.setParam('OutputFlag', False)
    model.update()
    return (model, variables)


def retrieveSolutionByIDs(solved_model, var_ids):
    """
    Retreive MILP solution by variable ids
    solver_model: gurobi solved model
    ids: aray of variable ids
    """
    x = []
    for var_id in ids:
        x.append(solved_model.getVarByName(var_id).x)
    return x


def isCandidatePair(X_sample, met_p, met_q):
    return (all(X_sample[met_q, :] <= X_sample[met_p, :])
            and not all(X_sample[met_q, :] == X_sample[met_p, :]))


def findCandidatePairs(model, n_samples=100):
    "Find candidate ordered concentration pairs"

    # Build sampling model
    e_constraints = {}
    e_vars = {'e_plus': {}, 'e_minus': {}}
    model.update()
    sampling_model = model.copy()
    sampling_model.setObjective(0, GRB.MINIMIZE)
    sampling_model.setParam('OutputFlag', False)

    # Get logx variables
    logx_vars = {}
    met_names = []
    for x in sampling_model.getVars():
        if 'logx' in x.VarName:
            logx_vars[x.VarName] = x
            met_names.append(x.VarName.replace('logx_', ''))

    for id in logx_vars.keys():

        # Add e_plus, e_minus variables
        var_id = f'e_plus_{id.replace("logx_", "")}'
        e_plus = sampling_model.addVar(lb=0, ub=1e6,
                                       obj=1.0, vtype=GRB.CONTINUOUS,
                                       name=var_id)

        var_id = f'e_minus_{id.replace("logx_", "")}'
        e_minus = sampling_model.addVar(lb=0, ub=1e6,
                                        obj=1.0, vtype=GRB.CONTINUOUS,
                                        name=var_id)
        e_vars['e_plus'][id] = e_plus
        e_vars['e_minus'][id] = e_minus

    sampling_model.update()

    # Add constraints
    for id in logx_vars.keys():
        c = f'e_vars["e_plus"]["{id}"] - e_vars["e_minus"]["{id}"] + logx_vars["{id}"] == 0'
        e_constraints[id] = sampling_model.addConstr(eval(c),
                                                     name=f'e_constr_{id}')

    # Sample solution space
    sampling_model.update()
    n_mets = len(logx_vars)
    lb, ub = np.log(par.x_min), np.log(par.x_max)
    X_sample = np.zeros((n_mets, n_samples))
    for n in range(n_samples):

        # print(f'sample number {n}')
        logx_rand_values = (ub - lb) * np.random.rand(n_mets) + lb

        for logx_rand, e_constr in zip(logx_rand_values, e_constraints.values()):
            e_constr.RHS = logx_rand

        sampling_model.update()
        sampling_model.optimize()
        X_sample[:, n] = retrieveSolutionByIDs(sampling_model, logx_vars.keys())

    candidatePairs = []
    for p in range(n_mets):
        for q in range(n_mets):
            if isCandidatePair(X_sample, p, q):
                candidatePairs.append((met_names[p], met_names[q]))

    return (X_sample, candidatePairs, sampling_model)
