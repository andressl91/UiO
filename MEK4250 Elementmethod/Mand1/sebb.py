from dolfin import *
import numpy as np
import matplotlib.pyplot as plt


def solver_1(N,k,l):
    mesh = UnitSquareMesh(N, N)

    V = FunctionSpace(mesh, "Lagrange", 2)

    class Right(SubDomain):
    	def inside(self, x, on_boundary):
    		return on_boundary and near(x[0],1)

    class Left(SubDomain):
    	def inside(self, x, on_boundary):
    		return on_boundary and near(x[0],0)


    right = Right()
    left = Left()
    bound = FacetFunction("size_t", mesh)
    bound.set_all(0)
    right.mark(bound, 1)
    left.mark(bound,2)
    #plot(bound); interactive()

    bc0 = DirichletBC(V, 0, right)
    bc1 = DirichletBC(V, 0, left)
    bcs = [bc0,bc1]

    u = TrialFunction(V)
    v = TestFunction(V)

    f = Expression("((pi*pi*k*k)+(pi*pi*l*l))*sin(pi*k*x[0])*cos(pi*l*x[1])",k=k,l=l)
    u_exact = interpolate(Expression("sin(pi*k*x[0])*cos(pi*l*x[1])",l=l,k=k),V)

    a = inner(grad(u), grad(v))*dx
    F = f*v*dx
    u_1 = Function(V)
    solve(a==F,u_1,bcs)
    return u_1,u_exact,V,mesh

    #print "mesh = %.d errornorm = %.8f  " %(N,errornorm(u_exact,u,norm_type = "H1" ,degree_rise = 3))
N = [2,4,8,16,32,64,128]
k = l = 100
b = np.zeros(len(N))
a = np.zeros(len(N))


for i in range(len(N)):
    u ,u_exact,V,mesh = solver_1(N[i],k,l)
    e = errornorm(u_exact,u,norm_type = "L2" ,degree_rise = 3)
    #u_ex = project(u_exact,V)
    u_norm = norm(u_exact,"H1")
    b[i] = np.log(e/u_norm)
    a[i] = np.log(1./(N[i]))
    #print u_norm
    #u_norm = sum(u_ex*u_ex)**0.5

A = np.vstack([a,np.ones(len(a))]).T
alpha_1,c_1 = np.linalg.lstsq(A,b)[0]
plt.plot(a,b,'o',label='Original data',markersize = 10)
plt.plot(a,alpha_1*a+c_1,'r',label='Fitted line')
plt.legend();plt.show()


B = np.zeros(2)
A = np.zeros((2,2))
A[0,0] = len(N)

for i in range(len(N)):
    A[0,1] += a[i]
    A[1,0] += a[i]
    A[1,1] += a[i]**2
    B[1] += a[i]*b[i]
    B[0] += b[i]
print A
c , alpha = np.linalg.solve(A, B)


print np.exp(c),alpha

plt.plot(a,b,'o',label='Original data',markersize = 10)
plt.plot(a,alpha*a+c,'r',label='Fitted line-my way')
plt.legend();plt.show()
#plt.plot(a,np.exp(alpha)*a+c); #plt.show()
#plt.plot(a,b); plt.show()

#print np.exp(x),y
