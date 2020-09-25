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


def buildLPModel(GEM):

    model = Model('ACR_model')
    GEM_rxn_ids = [rxn.id for rxn in GEM.reactions]
    S = cobra.util.array.create_stoichiometric_matrix(GEM, array_type='dense')

    # Add variables
    variables = {'logx': {}, 'v': {}, 'w': {}, 'logK': {}}
        
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
        
    for rxn in GEM.reactions:
        var_id = f'w_{rxn.id}'
#         x = model.addVar(lb=np.log(rxn.lower_bound), ub=np.log(rxn.upper_bound),
#                          obj=0.0, vtype=GRB.CONTINUOUS, name=var_id)
        x = model.addVar(lb=-GRB.INFINITY, ub=GRB.INFINITY,
                         obj=0.0, vtype=GRB.CONTINUOUS, name=var_id)
        variables['w'][var_id] = x
        
    for rxn in GEM.reactions:
        var_id = f'logK_{rxn.id}'
        x = model.addVar(lb=np.log(par.K_min), ub=np.log(par.K_max),
                         obj=0.0, vtype=GRB.CONTINUOUS, name=var_id)
        variables['logK'][var_id] = x

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

    for rxn in GEM.reactions:
        sum_str = ''
        for met in rxn.metabolites.keys():
            coeff = rxn.metabolites[met]
            if coeff < 0:
                sum_str += f'{-coeff} * variables["logx"]["logx_{met.id}"] +'
        sum_str = sum_str[:-1]
        if sum_str != '':
            c = f'variables["w"]["w_{rxn.id}"] == variables["logK"]["logK_{rxn.id}"] + {sum_str}'
            model.addConstr(eval(c), f'kinetics_{rxn_id}')

    # Adding log constraint w = log(v), only available in Gurobi 9
    for rxn in GEM.reactions:
        w = variables['w'][f'w_{rxn.id}']
        v = variables['v'][f'v_{rxn.id}']
        model.addGenConstrLog(v, w, name=f'w_{rxn.id} = log(v_{rxn.id})') 
#         model.addGenConstrExp(w, v, name=f'v_{rxn.id} = exp(w_{rxn.id})') 
    
    
    model.setParam('OutputFlag', par.gurobiFlag)
    model.update()
    return (model, variables)


def retrieveSolutionByIDs(solved_model, var_ids):
    """
    Retreive MILP solution by variable ids
    solver_model: gurobi solved model
    ids: aray of variable ids
    """
    x = []
    for var_id in var_ids:
        x.append(solved_model.getVarByName(var_id).x)
    return x


