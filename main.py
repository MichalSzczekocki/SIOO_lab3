import sys
from PyQt5.QtWidgets import QApplication, QWidget, QLabel, QVBoxLayout, QRadioButton, QComboBox, QSlider, \
    QPushButton
from PyQt5.QtCore import Qt
import numpy as np
import matplotlib.pyplot as plt
import warnings
import heapq

class Polynomial:

    def __init__(self, *coefficients):

        self.coefficients = list(coefficients)  # tuple is turned into a list

    def __repr__(self):

        return "Polynomial" + str(tuple(self.coefficients))

    def __str__(self):

        def x_expr(degree):
            if degree == 0:
                res = ""
            elif degree == 1:
                res = "x_1"
            else:
                res = "x_" + str(degree)
            return res

        degree = len(self.coefficients) - 1
        res = ""

        for i in range(0, degree + 1):
            coeff = self.coefficients[i]
            # nothing has to be done if coeff is 0:
            if abs(coeff) == 1 and i < degree:
                # 1 in front of x shouldn't occur, e.g. x instead of 1x
                # but we need the plus or minus sign:
                res += f"{'+' if coeff > 0 else '-'}{x_expr(degree - i)}"
            elif coeff != 0:
                res += f"{coeff:+g}{x_expr(degree - i)}"

        return res.lstrip('+')  # removing leading '+'

    def __call__(self, x):
        res = 0
        for coeff in self.coefficients:
            res = res * x + coeff
        return res

def repeatColumnNegative(Matrix, h):
    """Repeat column h multiplied by - 1"""
    (r, c) = Matrix.shape
    Matrix = np.hstack((Matrix[:, 0:h-1], -Matrix[:, [h-1]], Matrix[:, h-1:c]))
    return Matrix

def insertZeroToColumn(column, h):
    """insert zero to column"""
    k = np.size(column)
    column = np.vstack((column[0:h-1, [0]], np.array([[0]]), column[h-1:k, [0]]))
    return column

def minWithMask(x, mask):
    min = 0
    imin = 0
    n = np.size(x)
    for i in range(n):
        if mask[i] == 1:
            if min == 0:
                min = x[i]
                imin = i
            else:
                if min > x[i]:
                    min = x[i]
                    imin = i
    return (min, imin)

def SimplexSolution(type, A, B, C, D):
    """Calculates an optimal point for the linear programming model given by A*x <= B , Optimize z= C' * x
    Parametry funkcji:
    type -- type of optimization, it can be 'max' albo 'min'
    A    -- A matrix of the model (numpy array)
    B    -- B matrix of the model, column vector (numpy array)
    C    -- C matrix of the model, column vector (numpy array)
    D    -- column vector with the types of restrictions of the model (numpy array), 1 is <=, 0 is =, -1 is >=
            for <= restrictions do nothing
            for = add artificial variables
            for >= restrictions add slack variables and artificial variables in the auxilliary objective
            function (max --> -1 cost and 0 for other variables)
    """
    # m -- number of restrictions
    # n -- number of variables
    (m, n)= A.shape

    identityMatrix = np.eye(m)
    #Base = np.zeros(shape=(m))
    basic_vars = []
    count = n
    z_aux = [] #Auxilliary objective function
    a_vars= [] #Artificial variables
    vals = B
    BFS = [] #Dopuszczalne rozwiązanie bazowe

    for i in range(m):
        #identityMatrix[i,i] = identityMatrix[i,i]*D[i]
        if D[i]==1:
            #if the constraint is in the form of <=
            count = count + 1
            basic_vars = basic_vars + [count-1]
            a_vars = [a_vars, 0]

            if type == 'min':
                return
            else:
                return
        elif D[i]==0:
            #if the constraint is in the form of =
            #a_vars.append(1)
            count = count + 1
            basic_vars = basic_vars + [count-1]
            a_vars = [a_vars, 1]
            if type == 'min':
                return
            else:
                return
        elif D[i]==-1:
            #if the constraint is in the form of >=
            #a_vars.append(1)

            identityMatrix = repeatColumnNegative(identityMatrix, count + 1 - n)
            vals = insertZeroToColumn(vals, count + 1 - n)
            count = count + 2
            basic_vars = basic_vars + [count-1]
            a_vars = [a_vars, 0, 1]
            if type == 'min':
                return
            else:
                return
        else:
            print("Niewłaściwy przypadek")

    X = np.vstack((np.zeros((n, 1)), vals))
    A = np.hstack((A, identityMatrix))
    tableau = np.vstack((np.hstack((-np.transpose(C), np.array([[0]]))), np.hstack((A, B))))
    (rows, cols) = tableau.shape

    print('\nTablica Simpleks\n')
    print(tableau)
    print('\nZmienne bazowe\n')
    print(basic_vars)
    print('\nOptymalny punkt\n')
    print(X)

    iteration = 0
    
if __name__ == "__main__":
    #PRZYKLAD 5
    (z, x) = SimplexSolution('min', np.array([[3, 5], [5, 2]]),
                            np.array([[15], [10]]),
                            np.array([[-5], [-3]]),
                            np.array([[1], [1]]))