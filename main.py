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

def repeatColumnNegative(Matrix, h):
    (r, c) = Matrix.shape
    Matrix = np.hstack((Matrix[:, 0:h-1], -Matrix[:, [h-1]], Matrix[:, h-1:c]))
    return Matrix

def insertZeroToColumn(column, h):
    k = np.size(column)
    column = np.vstack((column[0:h-1, [0]], np.array([[0]]), column[h-1:k, [0]]))
    return column

def SimplexSolution(type, A, B, C, D):
    """FUNKCJA TO ZNALEZNIENIA OPTYMALNEGO ROZWIĄZANIA ZADANIA A*x <= B , Optymalizacja Z= C' * X
    Parametry funkcji:
    type -- type of optimization, it can be 'max' albo 'min'
    A    -- Macierz A współczynników ograniczeń
    B    -- Macierz B RHS ograniczeń modelu, wektor kolumnowy 
    C    -- Macierz C współczynników funkcji celu
    D    -- Wektor kolumnowy zawierający rodzaj znaków ograniczeń w zadaniu, 1 to <=, 0 to =, -1 to >=
    """
    # m -- ilość ograniczeń
    # n -- ilość zmiennych
    (m, n)= A.shape

    basic_vars = []
    count = n

    # MACIERZ Z NOWYMI ZMIENNYMI
    identityMatrix = np.eye(m)

    # WARTOŚCI NOWYCH ZMIENNYCH
    #Dopuszczalne rozwiązanie bazowe, basic feasible solution
    P = B

    # WSKAŹNIK POZYCJI ZMIENNYCH SZTUCZNYCH
    artificial= []

    #constraints = D #każdy kolejny element odpowiada równaniom ograniczającym, 1 to <=, 0 to =, -1 to >=

    for i in range(m):
        if D[i] == 1:

            C = np.vstack((C, [[0]]))
            count = count + 1
            basic_vars = basic_vars + [count-1]
            artificial = [artificial, 0]

        elif D[i] == 0:
            if type == 'min':
                C = np.vstack((C, [[100]]))
            else:
                C = np.vstack((C, [[-100]]))

            count = count + 1
            basic_vars = basic_vars + [count-1]

            artificial = [artificial, 1]
        elif D[i] == -1:
            # DODANIE ZMIENNYCH LUZU I SZCZTUCZNYCH
            if type == 'min':
                C = np.vstack((C, [[0], [100]]))
            else:
                C = np.vstack((C, [[0], [-100]]))

            identityMatrix = repeatColumnNegative(identityMatrix, count + 1 - n)
            P = insertZeroToColumn(P, count + 1 - n)

            # regist the artificial variable as basic variable
            count = count + 2
            basic_vars = basic_vars + [count-1]

            artificial = [artificial, 0, 1]
        else:
            print("Niewłaściwy przypadek dla ograniczeń")

    X = np.vstack((np.zeros((n, 1)), P))

    # DODANIE NOWEJ ZMIENNEJ 
    A = np.hstack((A, identityMatrix))

    # MACIERZ SIMPLEKS
    st = np.vstack((np.hstack((-np.transpose(C), np.array([[0]]))), np.hstack((A, B))))
    (rows, cols) = st.shape

    print('\nTablica Simpleks\n')
    print(st)
    print('\nZmienne bazowe\n')
    print(basic_vars)
    print('\nOptymalny punkt\n')
    print(X)

    # SPRAWDZENIE Z != 0 (GDY SĄ ZMIENNE BAZOWE)
    z_optimal = np.matmul(np.transpose(C), X)

    print('\nAktualna wartość Z funkcji celu\n\n', z_optimal)

    if z_optimal != 0:
        for i in range(m):
            if D[i] == 0 or D[i] == -1:
                if type == 'min':
                    st[0,:]= st[0,:] + 100 * st[1+i,:]
                else:
                    st[0,:]= st[0,:] - 100 * st[1+i,:]

        print('\nSkorygowana tablica simpleks\n')
        print(st)

    iteration = 0
    while True:
        if type == 'min':
            w = np.amax(st[0, 0:cols-1])
            iw = np.argmax(st[0, 0:cols-1])
        else:
            w = np.amin(st[0, 0:cols-1])
            iw = np.argmin(st[0, 0:cols-1])

        if w <= 0 and type == 'min':
            print('\nGlobalny punkt optimum\n')
            break
        elif w >= 0 and type == 'max':
            print('\nGlobalny punkt optimum\n')
            break
        else:
            iteration = iteration + 1

            print('\n----------------- Iteracja {} -------------------\n'.format(iteration))

            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                T = st[1:rows, cols-1] / st[1: rows, iw]

            R = np.logical_and(T != np.inf, T > 0)
            (k, ik) = minWithMask(T, R)

            # AKTUALNY ELEMENT Z WIERSZA
            cz = st[[0],:]

            # ELEMENT CENTRALNY
            pivot = st[ik+1, iw]

            # ELEMENT CENTRALNY PODZIELONY PRZEZ SAMEGO SIEBIE
            prow = st[ik+1,:] / pivot

            st = st - st[:, [iw]] * prow
            st[ik+1,:]= prow

            # NOWA ZMIENNA BAZOWA
            basic_vars[ik] = iw

            print('\nZmienne bazowe\n')
            print(basic_vars)

            basic = st[:, cols-1]
            X = np.zeros((count, 1))
            t = np.size(basic_vars)

            for k in range(t):
                X[basic_vars[k]] = basic[k+1]

            print('\nOptymalny punkt\n')
            print(X)
            C = -np.transpose(cz[[0], 0:count])
            z_optimal = cz[0, cols-1] + np.matmul(np.transpose(C), X)
            st[0, cols-1] = z_optimal

            print('\nTablica simpleks\n\n')
            print(st)

            print('\nAktualna wartość Z funkcji celu\n\n')
            print(z_optimal)

    # SPRAWDZENIE CZY KTÓRAŚ ZE ZMIENNYCH BAZOWYCH JEST DODATNIA
    tv = np.size(artificial)
    for i in range(tv):
        if artificial[i] == 1:
            if X[n + i] > 0:
                print('\nRozwiązanie nieoptymane\n')
                break

    return (z_optimal[0, 0], X)

if __name__ == '__main__':
    np.set_printoptions(suppress=True)
    #POKOLEI ODKOMENTOWYWAĆ PRZYKŁADY
    #PRZYKŁAD 2
    (z, x) = SimplexSolution('min', np.array([[1, 0, 0], [0, 1, 0], [0, 0, 5]]), np.array([[20], [6], [15]]), np.array([[5], [4], [6]]), np.array([[0], [0], [0]]))
    #PRZYKLAD 5
    #(z, x) = SimplexSolution('min', np.array([[3, 5], [5, 2]]), np.array([[15], [10]]), np.array([[-5], [-3]]), np.array([[1], [1]]))
    