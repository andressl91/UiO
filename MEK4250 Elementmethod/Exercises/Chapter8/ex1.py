#Estimiate ratio of non-zeros per unknown of the siffness matrix on the unit square
#considering u - laplace (u) = f, and stiffness matrix is A from Ax=b
from dolfin import *

Ratio = []
for i in [10]:
    for d in [1, 2, 3, 4]:
        mesh = UnitSquareMesh(i, i)
        mu  = 1

        V = FunctionSpace(mesh, 'Lagrange', d)

        u = TrialFunction(V)
        v = TestFunction(V)

        a = u*v*dx + inner(grad(u), grad(v))*dx
        A = assemble(a)

        #Number of unknowns
        unknowns = A.size(0)
        print '---------------------------------'
        print "Dimention %d" % d
        print "Matrix is a %d x %d matrix" % (unknowns, unknowns)
        print "With %d unknowns" % unknowns
         #Number of zeros
        zeros = A.nnz()
        print "Number of zeros in matrix %d" % zeros

        Ratio.append(zeros/unknowns)
print '---------------------------------'
print 'Ratio of non-zeros per unknown'
print Ratio
