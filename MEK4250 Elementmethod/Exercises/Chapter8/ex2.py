#Compute smallest and largest eigenvalues of the mass matrix and the
#stiffness matrix in 1D, 2D and 3D.
#Compute eigenvalues by from the eigenvalue problem
#Ax = lambda*M*x, where M is the mass matrix (FEM identity matrix)
from dolfin import *
import numpy as np

D = [1,2,3]
N = [10, 20]

for d in D:
    p = []; h = []
    r = np.zeros(len(N)-1)
    for n in N:
        if d == 1:
            mesh = IntervalMesh(n, 0, 1)
        elif d==2:
            mesh = UnitSquareMesh(n, n)
        else:
            mesh = UnitCubeMesh(n, n, n)


        V = FunctionSpace(mesh, 'Lagrange', 1)

        u = TrialFunction(V)
        v = TestFunction(V)

        a = u*v*dx + inner(grad(u), grad(v))*dx
        A = assemble(a)
        M = Identity(A.size(0))
        e = np.linalg.eigvals(A.array())
        e.sort()
        #VALUES
        p.append(e[-1]/e[0])
        h.append(mesh.hmin())

        #print e[0], e[-1]
    p = np.asarray(p); h = np.asarray(h)
    r[:] = np.ln(p[1:] - p[:-1])/ np.ln(h[1:] - h[:-1] )
    print r
