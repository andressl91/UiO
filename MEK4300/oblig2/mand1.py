from dolfin import *
from numpy import zeros
import matplotlib.pyplot as plt

L = 6
mesh = IntervalMesh(50, 0, L) 
V = FunctionSpace(mesh,'CG', 1)   
VV = V*V
vf, vh = TestFunctions(VV) #different testfunc for both equations

xeta = Function(V)
xeta = interpolate(Expression("x[0]"), V) 
xpoint  = xeta.vector().array()
#Remember that x[0] = eta....
xd = xpoint[1]-xpoint[0]

def newton(beta):

	bc0 = DirichletBC(VV, Constant((0,0)), "std::abs(x[0]) < 1e-10")
	bc1 = DirichletBC(VV.sub(1), Constant(1), "std::abs(x[0]-%d) < 1e-10" % L)
	#Sub 1 is the space for H, where f' = 1 at infinite

	#Guess first solution
	fh_ = interpolate(Expression(("x[0]", "1/L*x[0]"), L = L), VV)

	#L = L due to state what L is in Expression
	f_, h_ = split(fh_)

	F = h_*vf*dx -f_.dx(0)*vf*dx - inner(grad(h_),grad(vh))*dx + f_*h_.dx(0)*vh*dx + 		beta*(1-h_*h_)*vh*dx

	solve(F==0, fh_, [bc0,bc1])

	f_, h_ = fh_.split(True) #Splitting solutions for f and h 
	hp  = h_.vector().array()	#Solution of h derivative
	
	hder = zeros(len(xpoint)) #h derivative
	for i in range(len(hder)-1):
		hder[i] = (hp[i+1]- hp[i])/xd	
	
	return hp, hder

beta = [1, 0.3, 0, -0.1, -0.18,-0.19884] 

def plot():
	for i in range(len(beta)):
		hp, hder = newton(Constant(beta[i]))

		fig1 = plt.figure(1)
		ax1 = fig1.add_subplot(111)
		ax1.plot(xpoint, hp)

		fig2 = plt.figure(2)
		ax2 = fig2.add_subplot(111)
		ax2.plot(xpoint, hder)
	
	ax1.set_xlabel("$\eta$")
	ax1.set_ylabel("f derivative")	
	ax1.legend(["Beta = %f" % b for b in beta])	
	ax2.set_xlabel("$\eta$")
	ax2.set_ylabel("h derivative")	
	ax2.legend(["Beta = %f" % b for b in beta])

def negative():
	hp, hder = newton(Constant(beta[3]))
	fig3 = plt.figure(3)
	ax3 = fig3.add_subplot(111)
	ax3.plot(xpoint, hder)
	#Remember to change the bc1 to a negative value

#negative()
plot()
plt.show()	



