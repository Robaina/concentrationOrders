"""
The Gurobi python interface does not allow to build models in matrix form unlike both the
R and MATLAB interfaces. Instead constraints and variables have to be added iteratively.
This class allows to construct gurobi models in python from numpy arrays of constraints
which is oftentimes a more convenient formulation.
"""
import numpy as np
from gurobipy import Model, GRB

__author__ = "Semidan Robaina Estevez"


def updateObjective(model, cvec, sense='min'):

    if sense.lower() in 'min':
        obj_sense = GRB.MINIMIZE
    else:
        obj_sense = GRB.MAXIMIZE

    x = model.getVars()
    o = ''
    for i, coeff in enumerate(cvec):
        if coeff != 0:
            o += 'x[' + str(i) + '] * ' + str(coeff) + '+'
    o = o[:-1]

    model.setObjective(eval(o), obj_sense)
    model.update()
    return model


def updateRHS(model, rhs):
    """
    This version update the last set of constraints with size equal to RHS
    """
    N = len(rhs)
    Constrs = model.getConstrs()[-N:]
    for i, constr in enumerate(Constrs):
        constr.setAttr('rhs', rhs[i])

    model.update()
    return model


def constructGurobiModel(c, G, h, A, b, lb=None, ub=None, sense='min',
                         binaryVariables=None, variableNames=None, modelName=None):
    """
    Construct gurobipy model from a matrix formulation

    Arguments:
    ---------
    c: 1D array, the objective vector
    A: 2D array, the constraint matrix
    lb, ub: 1D array, lower and upper bounds for the variables
    sense: str, the optimization sense, 'min' or 'max'
    binaryVariables: array (optional), column indices of A corresponding
                     to binary variables
    variableNames: list of str (optional), the names of the variables
    """
    N_eq_constraints, N_vars = np.shape(A)
    # N_ineq_constraints = len(h)

    if modelName is None:
        model_name = 'model'
    else:
        model_name = modelName

    if variableNames is None:
        var_names = ['x' + str(n) for n in range(N_vars)]
    else:
        var_names = variableNames

    # Create variable types
    var_types = [GRB.CONTINUOUS for v in range(N_vars)]
    if binaryVariables is not None:
        for i in binaryVariables:
            var_types[i] = GRB.BINARY

    # Create model object and add variables
    m = Model(model_name)
    x = m.addVars(range(N_vars), lb=-GRB.INFINITY, ub=GRB.INFINITY,
                  obj=c, vtype=var_types, name=var_names)

    # Add inequality constraints
    for i, row in enumerate(G):
        s = ''
        for j, coeff in enumerate(row):
            if coeff != 0:
                s += 'x[' + str(j) + '] * ' + str(coeff) + '+'
        s = s[:-1]
        s += ' <= ' + str(h[i])
        m.addConstr(eval(s))

    # Add equality constraints
    for i, row in enumerate(A):
        s = ''
        for j, coeff in enumerate(row):
            if coeff != 0:
                s += 'x[' + str(j) + '] * ' + str(coeff) + '+'
        s = s[:-1]
        s += ' == ' + str(b[i])
        m.addConstr(eval(s))

    # Set objective
    m.update()
    m = updateObjective(m, c, sense='min')

    return m
