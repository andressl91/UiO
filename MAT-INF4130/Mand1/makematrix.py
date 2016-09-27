
import numpy as np
from scipy.linalg import hilbert

def hilb(i, j):
    return 1 / ((i+1) + (j+1) - 1)


n = 3
Pnum = np.eye(n+1)

Hilbert = np.fromfunction(lambda i, j: 1 / ((i+1) + (j+1) - 1), (n+1, n+1), dtype=float)

altsign = [(-1)**i for i in range(n+1)]
altsign = np.asarray(altsign)
altsignmatrix = np.matrix([altsign*(-1)**i for i in range(n+1)])
Nnum = Hilbert + np.multiply(altsignmatrix, Hilbert)
print Nnum
#print Nnum
#for i in range(1, n):
#    print "i = %d" % i
#    for j in range(i):
#        print Nnum[j, j]
        #print  (-1*(Nnum[j, i]))
        #Pnum[j, i] = Nnum[j, j]/(-1*(Nnum[j, i]))
#for i in range(n+1):
#print Pnum

for i in range(1, n+1):
    print "I is %d" % i
    #print Nnum[:i, :i]
    #print  (-1*(Nnum[:i, i]))
    print np.divide( Nnum[:i, :i],(-1*(Nnum[:i, i])))


#Symbolic
from sympy import *
from mpmath import hilbert
#from sympy.matrices import *
from IPython.display import display

def f(i, j):
    if i % 2 == 0:
        return 1 if j % 2 == 0 else -1
    if i % 2 == 1:
        return -1 if j % 2 == 0 else 1

altsignmatrix = Matrix(n+1, n+1, f)
mat = Matrix(n+1, n+1, hilb)
m = eye(n+1)
for i in range(n+1):
    for j in range(n+1):
        m[i, j] = mat[i, j] + altsignmatrix[i,j]*mat[i, j]
#display(m)
