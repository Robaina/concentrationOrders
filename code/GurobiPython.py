"""
The Gurobi python interface does not allow to build models in matrix form unlike both the
R and MATLAB interfaces. Instead constraints and variables have to be added iteratively.
This class allows to construct gurobi models in python from numpy arrays of constraints
which is oftentimes a more convenient formulation.

Semidan Robaina Estevez
"""
import numpy as np
from gurobipy import Model, GRB

class GurobiModel:
    """
    Constructs a gurobipy model from a matrix formulation. Currently only LPs
    and MILPs are covered, i.e., the general optimization problem:
                  min   cx
                  s.t.
                      Ax <= b (>=, ==)
                      lb <= x <= ub
    Arguments:
    ---------
    c: 1D array, the objective vector
    A: 2D array, the constraint matrix
    lb, ub: 1D array, lower and upper bounds for the variables (default, 0 to Inf)
    modelSense: str, the optimization sense, 'min' or 'max'
    sense: array-like of str, the constraints sense: <=, == or >= (default <=)
    binaryVariables: array (optional), column indices of A corresponding
                     to binary variables (default continuous)
    variableNames: array-like of str (optional), the names of the variables
    modelName: str, (optional) the name of the gurobi model.
    """

    def __init__(self, c, A, b, lb=0, ub=GRB.INFINITY, modelSense='min',
                 sense=None, binaryVariables=None,
                 variableNames=None, modelName='model'):
        self.obj = c
        self.A = A
        self.rhs = b
        self.ub = ub
        self.lb = lb
        self.modelSense = modelSense
        self.binaryVariables = binaryVariables
        if modelName is None:
            self.modelName = 'model'
        else:
            self.modelName = modelName

        self.nConst, self.nVars = np.shape(A)
        if sense is None:
            self.sense = ['<=' for _ in range(self.nConst)]
        else:
            self.sense = sense
        self.varType = np.array(['C' for _ in range(self.nVars)])
        self.varType[self.binaryVariables] = 'B'
        if variableNames is None:
            self.varNames = ['x' + str(n) for n in range(self.nVars)]
        else:
            self.varNames = variableNames

    def construct(self):
        """
        Builds a gurobipy model object
        """
        model = Model(self.modelName)

        x = model.addVars(range(self.nVars), lb=self.lb, ub=self.ub,
                          obj=self.obj, vtype=self.varType, name=self.varNames)

        for i, row in enumerate(self.A):
            s = ''
            for j, coeff in enumerate(row):
                if coeff != 0:
                    s += 'x[' + str(j) + '] * ' + str(coeff) + '+'
            s = s[:-1]
            s += ' ' + self.sense[i] + ' ' + str(self.rhs[i])
            model.addConstr(eval(s))
        
        model = self.updateObjective(model, self.obj)
        return model

    @staticmethod
    def updateObjective(model, c, sense='min'):
        """
        Updates the objective vector c of a linear model
        """
        if sense.lower() in 'min':
            objsense = GRB.MINIMIZE
        else:
            objsense = GRB.MAXIMIZE
        model.update()
        
        x = model.getVars()
        o = ''
        for i, coeff in enumerate(c):
            if coeff != 0:
                o += 'x[' + str(i) + '] * ' + str(coeff) + '+'
        o = o[:-1]

        model.setObjective(eval(o), objsense)
        model.update()
        return model
        
    @staticmethod
    def updateRHS(model, b):
        """
        Updates the right-hand-side vector b of the model constraints
        """
        model.update()
        Constrs = model.getConstrs()
        for i, constr in enumerate(Constrs):
            constr.setAttr('rhs', b[i])

        model.update()
        return model