#Waveproject
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter
from mpl_toolkits.mplot3d import Axes3D

import numpy as np
import sympy as sym
import matplotlib.pyplot as plt
import sys, unittest

def symbolicPDE(u_e,qval):
	x,y,t,q,b,A,B,w,kx,ky= sym.symbols('x,y,t,q,b,A,B,w,kx,ky')  # global symbols
	u = u_e; q = qval
	u_x = sym.diff(u,x)
	u_y = sym.diff(u,y)
	q_x = sym.diff(q,x)
	q_y = sym.diff(q,y)
	rhs = q_x*u_x + q*sym.diff(u,x,x) + q_y*u_y + q*sym.diff(u,y,y)
	lhs = sym.diff(u,t,t) + b*sym.diff(u,t)
	return lhs-rhs

def init(x,y,t, test):
	if test == "constant":
		return 2.0
	if test == "plug":
		return 1.0
	else:
		return np.cos(np.pi/1.0*x)*np.cos(np.pi/1.0*y)

def solver(version, dt, dx, dy, Lx, Ly, T, b, V, I, qfunc, exact, f, plot):
	 # fixe
	beta = 0.9
	Nx = int(round(Lx/dx))
	Ny = int(round(Ly/dy))
	x = np.linspace(0, Lx, Nx+1)
	y = np.linspace(0, Ly, Ny+1)
	xv = x[:,np.newaxis]
	yv = y[np.newaxis,:]
 	y = np.linspace(0, Ly, Ny+1)
 	X, Y = np.meshgrid(x, y)

	Nt = int(round(T/float(dt)))
	t = np.linspace(0, Nt*dt, Nt+1) # mesh points in time
	errors = np.zeros(len(t)) #For storing errors in the scheme
	
	u = np.zeros((Nx+3,Ny+3))
	u_1 = np.zeros((Nx+3,Ny+3))
	u_2 = np.zeros((Nx+3,Ny+3))
	q = np.zeros((Nx+3,Ny+3))
	to_return = np.zeros((Nx+1,Ny+1))
	
	q[1:-1,1:-1] = qfunc(xv[:],yv[:])
	q[0,1:-1] = q[2,1:-1]; q[-1,1:-1] = q[-3,1:-1]
	q[1:-1,0] = q[1:-1,2]; q[1:-1,-1] = q[1:-1,-3]

	
	stability_limit = (1/float(np.amax(q)))*(1/np.sqrt(1/dx**2 + 1/dy**2))
	if dt <= 0: # max time step?
		safety_factor = -dt # use negative dt as safety factor
		dt = safety_factor*stability_limit
	elif dt > stability_limit:
		print 'error: dt=%g exceeds the stability limit %g' % \
				(dt, stability_limit)
	
	#Initial value
	
	u_2[1:-1,1:-1] = I(xv,yv)
	u_2[0,1:-1] = u_2[2,1:-1]; u_2[-1,1:-1] = u_2[-3,1:-1]
	u_2[1:-1,0] = u_2[1:-1,2]; u_2[1:-1,-1] = u_2[1:-1,-3]

	errors[0] = np.max(abs(exact(xv[:],yv[:],t[0])-u_2[1:-1,1:-1]))

	#First value
	u_1[1:-1,1:-1] = 0.5*(2*u_2[1:-1,1:-1]+2*dt*V*(dt*b/2.0-1)+ \
		(0.5*(dt/dx)**2*(
		(q[2:,1:-1]+q[1:-1,1:-1])*(u_2[2:,1:-1]-u_2[1:-1,1:-1])-(q[1:-1,1:-1]+q[:-2,1:-1])*(u_2[1:-1,1:-1]-u_2[:-2,1:-1])) +
		0.5*(dt/dx)**2*(
		(q[1:-1,2:]+q[1:-1,1:-1])*(u_2[1:-1,2:]-u_2[1:-1,1:-1])-(q[1:-1,1:-1]+q[1:-1,:-2])*(u_2[1:-1,1:-1]-u_2[1:-1,:-2])))
		+dt**2*f(xv[:],yv[:],t[1]))

	u_1[0,1:-1] = u_1[2,1:-1]; u_1[-1,1:-1] = u_1[-3,1:-1]
	u_1[1:-1,0] = u_1[1:-1,2]; u_1[1:-1,-1] = u_1[1:-1,-3]

	errors[1] = np.max(abs(exact(xv[:],yv[:],t[1])-u_1[1:-1,1:-1]))
	if version == "scalar":
		for k in range(2,Nt):
			for i in range(Nx):
				for j in range(Ny):
					u[i][j] = 1.0/(1+0.5*b*dt)*(2*u_1[i,j]+u_2[i,j]*(0.5*dt*b-1) \
					+ (0.5*(dt/dx)**2*(
					(q[i+1:,j]+q[i,j])*(u_1[i+1:,j]-u_1[i,j])-(q[i,j]+q[i-1,j])*(u_1[i,j]-u_1[i-1,j])) +
					0.5*(dt/dx)**2*(
					(q[i,j+1]+q[i,j])*(u_1[i,j+1]-u_1[i,j])-(q[i,j]+q[i,j-1])*(u_1[i,j]-u_1[i,j-1])))
					+dt**2*f(x[i],y[j],t[k]))

			u_2[:] = u_1[:]
			u_1[:] = u[:]

			if plot == True:
				num_frames = 80
				skip_frame = int(len(t)/float(num_frames))
				if i % skip_frame == 0 or i == t.size-1:
					fig = plt.figure()
					fig.suptitle('Plot' , fontsize=14, fontweight='bold')
					ax = fig.add_subplot(111, projection='3d')
					ax.set_title("Time %.1f seconds" % (i*dt))
	
					ax.set_xlim([0, Lx]); ax.set_ylim([0, Ly])
					ax.set_zlim([0, 3])	
					ax.plot_wireframe(X, Y, u_2[1:-1,1:-1], rstride=2, cstride=2, color="blue")
					#ax.plot_wireframe(X, Y,exact(xv[:],yv[:],t[i-1]) , rstride=2, cstride=2,color="g")
					ax.view_init(15,0)
					ax.set_xlabel('x values')
					ax.set_ylabel('y values')
					plt.savefig('frame_%04d.png' % i)
					plt.close(fig)
		
	elif version == "vectorized":
		for i in range(2,Nt):
			u[1:-1,1:-1] = 1.0/(1+0.5*b*dt)*(2*u_1[1:-1,1:-1]+u_2[1:-1,1:-1]*(0.5*dt*b-1) \
			+ (0.5*(dt/dx)**2*(
			(q[2:,1:-1]+q[1:-1,1:-1])*(u_1[2:,1:-1]-u_1[1:-1,1:-1])-(q[1:-1,1:-1]+q[:-2,1:-1])*(u_1[1:-1,1:-1]-u_1[:-2,1:-1])) +
			0.5*(dt/dx)**2*(
			(q[1:-1,2:]+q[1:-1,1:-1])*(u_1[1:-1,2:]-u_1[1:-1,1:-1])-(q[1:-1,1:-1]+q[1:-1,:-2])*(u_1[1:-1,1:-1]-u_1[1:-1,:-2])))
			+dt**2*f(xv[:],yv[:],t[i]))

			u[0,1:-1] = u[2,1:-1]; u[-1,1:-1] = u[-3,1:-1]
			u[1:-1,0] = u[1:-1,2]; u[1:-1,-1] = u[1:-1,-3]
		
			u_2[:] = u_1[:]
			u_1[:] = u[:]
			#Tenk riktig pa timesteppet fordi vi starter paa 2, saa losning er for 3
			errors[i] = np.max(np.abs(exact(xv[:],yv[:],t[i])-u_1[1:-1,1:-1]))

			if plot == True:
				num_frames = 90
				skip_frame = int(len(t)/float(num_frames))
				
				if i % skip_frame == 0 or i == t.size-1:
					fig = plt.figure()
					fig.suptitle('Wave plot' , fontsize=14, fontweight='bold')
					ax = fig.add_subplot(111, projection='3d')
					ax.set_title("Time %.1f seconds" % (i*dt))
					
					ax.set_xlim([0, Lx]); ax.set_ylim([0, Ly])
					ax.set_zlim([0, 3])		
					ax.plot_wireframe(X, Y, u_2[1:-1,1:-1], rstride=2, cstride=2, color="blue")
					#ax.plot_wireframe(X, Y,exact(xv[:],yv[:],t[i-1]) , rstride=2, cstride=2,color="g")
					ax.view_init(25,0)
					ax.set_xlabel('x values')
					ax.set_ylabel('y values')
					plt.savefig('frame_%04d.png' % i)
					plt.close(fig)
		bigerr = np.max(errors)
		timeloc = np.argmax(errors)
		return bigerr 
		

def convergence_rates(version, solver_function, ue, m, h,Lx, Ly, T, b, V, I, qfunc, f,plot=False):
		     			
    dt_values = []
    E_values = []
    for i in range(m):
    	dt = h/10.
    	Nx = 1.0/h
    	Ny = 1.0/h
    	dx = float(Lx/(Nx))
    	dy = float(Ly/(Ny))

       	u_quad_diff = solver_function(version, dt, dx, dy, Lx, Ly, T, b, V, I, qfunc, ue, f, plot)
        #E = np.sqrt(dt*u_quad_diff)
        E = u_quad_diff
        dt_values.append(dt)
        E_values.append(E)
        h = h/2.0
        
      
    r = [np.log(E_values[i-1]/E_values[i])/
         np.log(dt_values[i-1]/dt_values[i])
         for i in range(1, m, 1)]
    return r, dt_values, E_values
	

class main():
	#f=0,V=0, b and q are free, I = u = c
	def __init__(self):
		self.dx = None; self.dy = None
		self.dt = None; self.q = None
		self.I = None; 
		self.f = None

	def const(self, version, plot):
		dx = 0.05; dy = 0.05; self.dx = dx; self.dy = dy
		b = 3; V = 0.0
		Lx=1; Ly=1; T=4
		I = np.vectorize(lambda x,y: 1)
		exact = np.vectorize(lambda x,y,t: 1)
		q = np.vectorize(lambda x,y: 1)
		f = np.vectorize(lambda x,y,t: 0)
		self.I = I; self.q = q; self.f = f

		c = np.sqrt(q(1,1))
		dt = (1/float(c))*(1/np.sqrt(1/dx**2 + 1/dy**2))
		return solver(version, dt, dx, dy, Lx, Ly, T, b, V, I, q, exact, f, plot)
	#Exercise 3.3
	def plug(self,version,plot):
		dx = 0.03; dy = 0.03
		self.dx = dx; self.dy = dy; 
		b = 0; V = 0.0
		Lx=4; Ly=4; T=4
		test = "plug"; plot = True
		I = np.vectorize(lambda x,y: 0 if abs(x-Lx/2.0) > 0.08 else 1)
		q = np.vectorize(lambda x,y: 2); self.q = q
		f = np.vectorize(lambda x,y,t: 0)
		exact = np.vectorize(lambda x,y,t: 0)
		c = np.sqrt(q(1,1))
		dt = dx/c
		self.dt = dt
		solver(version, dt, dx, dy, Lx, Ly, T, b, V, I, q, exact, f, plot)
	#Exercise 3.4
	def physprob(self, version, plot):
		dx = 0.05; dy = 0.05
		self.dx = dx; self.dy = dy
		b = 0.0; V = 0.0
		Lx=6.; Ly=6.; T=6
		
		plot = True
		A = 1; B = 1; mx = 1.0; my = 1.0
		kx = mx*np.pi/Lx; ky = my*np.pi/Ly
		#Im refelection loc, Is is witdth
		
		Nx = int(round(Lx/dx)); Ny = Nx
		x = np.linspace(0, Lx, Nx+1);
		y = np.linspace(0, Lx, Nx+1);
		xv = x[:,np.newaxis]
		yv = x[np.newaxis,:]
		X, Y = np.meshgrid(x, x)
		Im = 0; Is = 1.2; Io = 1; Ia = 1
		I = np.vectorize(lambda x,y: Io+Ia*np.exp(-((y-Im)/Is)**2))
		task = "box"

		if task == "dump":
			Bmx = Lx/3.; Bmy = Lx/3.; Bs = 1; Bo = -0.7; Ba = 1; scale = 8 	
			B = np.vectorize(lambda x,y,t: Bo+Ba*np.exp(-((y-Bmx)/Bs)**2 -((0-Bmy)/(scale*Bs)))**2)
		if task == "cos":
			Bmx = 1.; Bmy = 1.; Bs =1.; Bo = -0.17; Ba = 1.1; scale = 8
			B = np.vectorize(lambda x,y,t: Bo+Ba*np.cos(np.pi*(x-Bmx)/(2.*Bs))*np.cos(np.pi*(y-Bmy)/(2.*Bs)) \
								if 0 < np.sqrt((x-Lx/2.)**2+(y-Lx/2.)**2) <= Bs else Bo)
		if task == "box":
			Bmx = Lx/2.; Bmy = Lx/2.3; Bs = 0.7; Bo = -0.7; Ba = 1.6; scale = 1 	
			B = np.vectorize(lambda x,y,t: Bo+Ba if Bmx-Bs*2 <= x <= Bmx+Bs*2 and\
							 Bmy-scale*Bs <= y <= Bmy+scale*Bs else Bo)
		
		q = np.vectorize(lambda x,y: 9.81*(I(-1,7)-B(x,y,1)))
		f = np.vectorize(lambda x,y,t: 0)
		c = np.max(np.sqrt(q(x,y)))
		dt = (1/float(c))*(1.0/np.sqrt(1.0/dx**2 + 1.0/dy**2)) 
		
		solver(version, dt, dx, dy, Lx, Ly, T, b, V, I, q, B, f, plot)

	def standing(self, version, plot):
		dx = 0.1; dy = 0.1
		self.dx = dx; self.dy = dy
		b = 0.0; V = 0.0
		Lx=5; Ly=5; T=6
		
		test = "plug";
		A = 1; B = 1; mx = 2.0; my = 2.0
		kx = mx*np.pi/Lx; ky = my*np.pi/Ly

		I = np.vectorize(lambda x,y: A*np.cos(kx*x)*np.cos(ky*y) )
		q = np.vectorize(lambda x,y: 2)
		f = np.vectorize(lambda x,y,t: 0)
		c = np.sqrt(q(1,1))

		w = np.sqrt(q(1,1)*(kx**2+ky**2))
		ue = np.vectorize(lambda x,y,t: A*np.cos(kx*x)*np.cos(ky*y)*np.cos(w*t))

		m = 6; h = 1.0
		err = np.zeros(m+1)
		R,T,V = convergence_rates(version, solver, ue, m, h,Lx, Ly, T, b, V, I, q, f,plot=False)
		for i in range(len(R)):
			print "For dt = %.4f, the convergence rate is %.4f" % (T[i],R[i])
	
	def mansol(self, version, plot):

		def makef(Aval, Bval, wval, kxval, kyval, damping, qval):
			#Next session, implement w in the symbolicPDE mehtod
			x,y,t,A,B,w,kx,ky,q,b= sym.symbols('x,y,t,A,B,w,kx,ky,q,b')  # global symbols
			q = qval
			ue = (A*sym.cos(w*t)+B*sym.sin(w*t))*sym.exp(-sym.sqrt(q)*t)*sym.cos(kx*x)*sym.cos(ky*y)
			f = symbolicPDE(ue, q)
			f = f.subs([(A,Aval),(B,Bval),(kx,kxval),(ky,kyval),(b,damping)])
			w = sym.solve(f,w) #Takes to long to calculate
			w = w.subs([(A,Aval),(B,Bval),(kx,kxval),(ky,kyval),(b,damping)])
			f = np.vectorize(sym.lambdify((x,y,t), f))
			w = np.vectorize(sym.lambdify())
			return f,w

		dx = 0.05; dy = 0.05
		self.dx = dx; self.dy = dy
		V = 0.0; plot = False
		Lx=1; Ly=1; T=6
		b = 1 #Damping
		A = 1; B = 1; mx = 1.0; my = 1.0
		kx = mx*np.pi/Lx; ky = my*np.pi/Ly

		I = np.vectorize(lambda x,y: A*np.cos(kx*x)*np.cos(ky*y) )
		q = np.vectorize(lambda x,y: 1+x)
		#w = np.sqrt(q(1,1)*(kx**2+ky**2))
		c = np.sqrt(q(1,1))	; w = np.sqrt(q(1,1)*(kx**2+ky**2))
		ue = np.vectorize(lambda x,y,t: A*np.cos(kx*x)*np.cos(ky*y)*np.cos(w*t))
		x,y = sym.symbols('x,y')
		f,w = makef(A,B,w,kx,ky,b,x+1)
		dt = (1/float(c))*(1.0/np.sqrt(1.0/dx**2 + 1.0/dy**2)) 

		#print solver(version, dt, dx, dy, Lx, Ly, T, b, V, I, q, ue, f, plot)


class TestMySolver(unittest.TestCase):

	def setUp(self):
		self.prog = main()
		self.version = "vectorized"
		self.plot = False

	def test_constant_solution(self):
		x,y,t= sym.symbols('x,y,t'); constant = 1
		f = sym.lambdify((x,y,t), symbolicPDE(constant,x*y))
		randval=np.random.randint(10,size=3)
	   	x = randval[0]; y = randval[1]; t = randval[2];
	   	error = int(self.prog.const("vectorized",False)) 
	   	#Check f = 0
		self.assertEqual(0,f(x,y,t)) 
		#Check Numerical solution stays constant
		self.assertEqual(0,error) 
		#Check I = c
		self.assertEqual(constant,self.prog.I(x,y)) 

	def test_plug_solution(self):
		print "Will raise stability_limit error"
		self.prog.plug("vectorized",False)
		dt = self.prog.dt
		dx = self.prog.dx
		dy = self.prog.dy
		q = self.prog.q
		condx = np.sqrt(q(1,1))*dt/dx
		condy = np.sqrt(q(1,1))*dt/dy
		#Check wave velocitycondition
		self.assertEqual(1,condx)
		self.assertEqual(1,condy) 

if __name__ == '__main__':
	version = "vectorized"; plot = True
	#Setting plot to True, stores plots in 
	#your current folder. #Check of the class Main
	#run a method if you want. All results are presented
	#in the report. Here is a testrun for convergencerate
	#Exerice 3.4 
	run = main()
	run.standing(version, plot)
	#Unittest checks my conditions for the exercises stated in the reptort
	#Unittest needed
	#unittest.main()
	