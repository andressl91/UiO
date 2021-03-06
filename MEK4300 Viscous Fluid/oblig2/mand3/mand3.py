from numpy import arctan
from dolfin import *
import matplotlib.pyplot as plt

v_star = 0.05
kappa = 0.41
Re = 1000
A = 26
ht =  1 #Top wall
hb = -1 #Bottom wall
nu = v_star/Re 

mesh = IntervalMesh(250,0,1)
x = mesh.coordinates()
x[:,0] = ht -arctan(pi*(x[:,0])) / arctan(pi)  

V = FunctionSpace(mesh, 'CG', 1)
u = TrialFunction(V)
v = TestFunction(V)

#y = Expression("x[0]")
#y = interpolate(y, V)
#plot(y) For understanding
#interactive()

bcs = DirichletBC(V, 0, "std::abs(x[0]) < DOLFIN_EPS")

l = Expression("kappa*x[0]*(1- exp(-x[0]*v_star/(nu*A)))", kappa = kappa, nu = nu, A = A, v_star = v_star) 

F = -nu*inner(grad(u), grad(v))*dx + (v_star)**2/ht*v*dx
a = lhs(F)
L = rhs(F)

u_ = Function(V)
solve(a==L, u_, bcs)

F_ = nu*inner(grad(u_), grad(v))*dx - (v_star)**2/ht*v*dx + l*l*abs(u_.dx(0))*u_.dx(0)*v.dx(0)*dx

#Newtons method
solve(F_== 0, u_, bcs)
plot(u_)
interactive()

#Testing the theoretical values
y_star = x*Re
u_v = u_.compute_vertex_values()
u_star = u_v/v_star

plt.figure(1)
plt.plot(y_star, u_star, label = "$y^{+} = u^{+}$")
plt.axis([0, 5, 0, 7])
plt.title("Plotting $y^{+}$ against $u^{+}$ for $y^{+}$ < 5") 
plt.legend(loc="upper left")


B = 5.5
ut = interpolate(Expression("1/kappa * log(x[0]*Re) + B", B = B, kappa = kappa, Re = Re), V)
u_val = ut.compute_vertex_values()

plt.figure(2)
plt.plot(y_star, u_val,'b', label = "Theoretical")
plt.plot(y_star, u_star, 'g', label = "Numerical, B = %1.1f" % B)
plt.axis([30 ,1000, 0, 25])
plt.title("Comparing theoretical result for numerical result for $y^{+} > 30$") 
plt.legend(loc='upper left')

plt.show()

u_deri = project(u_.dx(0),V)
li = interpolate(l, V)
li = li.compute_vertex_values()


vt = li[0]**2*u_deri(1) #blir null, men fysikken sier at det er turb, spred, der, virvler
print vt











