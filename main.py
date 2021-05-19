from sys import stdin
import itertools
EPS = 1e-4

class Position:
    def __init__(self, column, row):
        self.column = column
        self.row = row

# Takes an int "n" that defines how many total inequalty equations there are
# Takes an int "m" that defines how many total variables there are
# Takes a list of lists "a" of size n x m that contains the inequalities.
# Takes a list "b" of size n of the maximums of each inequality
# Takes a list "c" that represents the optimization function to maximize
# For example:
# x + y - 3z <= 10
# -5x + 10y <= 50
# 3x - 2y -4z <= 9
# Maximize: x + 6y -3z
# n, m = 3, 3
# a = [[1,1,-3],[-5,10,0],[3,-2,-4]]
# b = [10, 50, 9]
# c = [1, 6, -3]
def ReadEquation():
    n, m = map(int, input().split())
    a = []
    for row in range(n):
        a.append(list(map(float, input().split())))
    b = list(map(float, input().split()))
    c = list(map(float, input().split()))
    return a, b, c, n, m

if __name__ == "__main__":

    exit(0)