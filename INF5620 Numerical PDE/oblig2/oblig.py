import numpy as np
import matplotlib.pyplot as plt
import sympy as sym
import os, glob
 #FIX DX TIL EXAMEN; SJEKK s 46
#os.remove('fig_*.png')
#convert -delay 10 -loop 0 frame_*.png animation.gif
def fvalues(task):
	x, t, L, w= sym.symbols('x, t, L, w')  # global symbols
	if task == "a":
		q = 1 + (x-L/2)**4
	else:
		q = 1 + sym.cos(sym.pi*x/L)

	u = sym.cos(x*sym.pi/L)*sym.cos(w*t)
	u_x = sym.diff(u,x)
	q_x = sym.diff(q,x)
	f = sym.diff(u,t,t) -q_x*u_x - q*sym.diff(u,x,x)
	sym.simplify(f) 
	sym.simplify(q)
	f = sym.lambdify((x,t,L,w), f)
	q = sym.lambdify((x,L), q)
	return f, q


def fv(x,t,L,w):
	return (-16*L**2*w**2*np.cos(np.pi*x/L) - 8*np.pi*L*(L - 2*x)**3*np.sin(np.pi*x/L) + \
	np.pi**2.0*((L - 2*x)**4 + 16)*np.cos(np.pi*x/L))*np.cos(t*w)/(16*L**2)

def qval(x, L):
	return 1.0 + (x-L/2.0)**4

def exact(x, t, L, w):
	return np.cos(np.pi*x/L)*np.cos(w*t)

def solver(V, w, dx, L, dt, T, task):
	Nt = int(round(T/dt))
	t = np.linspace(0, T, Nt+1)
	Nx = int(round(L/dx))
	x = np.linspace(0, L, Nx+1)

	c = dt/dx
	u = np.zeros(Nx+1)
	u_1 = np.zeros(Nx+1)#last
	u1  = np.zeros(Nx+1)#next
	
	
	fval, q = fvalues(task) #Test of sympy implementation, remove if needed
	fval = np.vectorize(fval)
	q = np.vectorize(q); q = q(x, L)

	#Initial Condition
	u_1[:] = exact(x, 0, L, w)
	#First step
	u[1:-1] = dt*V + u_1[1:-1]+c**2*0.25*\
	((q[1:-1]+q[2:])*(u_1[2:]-u_1[1:-1])-(q[1:-1]+q[:-2])*(u_1[1:-1]-u_1[:-2])) +\
	0.5*dt**2*fval(x[1:-1],0,L,w)

	if task == "c":
		u[0] = u[1]
		u[-1] = u[-2]
	if task == "d": 
		u[0] = dt*V + u_1[0]+c**2*0.25*(-(q[1]+q[0])*(u_1[1]-u_1[0])) + 0.5*fval(x[0],0,L,w)
		u[-1] = dt*V + u_1[-1]+c**2*0.25*(-(q[-1]+q[-2])*(u_1[-1]-u_1[-2])) + 0.5*fval(x[-1],0,L,w)
	else: #Check brok
		u[0] = dt*V + u_1[0]+c**2 *0.5*(q[0]+q[1])*(u_1[1] -u_1[0]) + 0.5*dt**2*fval(x[0],0,L,w)
		u[-1] = dt*V + u_1[-1]+c**2 *0.5*(q[-1]+q[-2])*(u_1[-2] -u_1[-1]) + 0.5*dt**2*fval(x[-1],0,L,w)

	error = u
	value = np.max(abs(u-exact(x, t[1], L, w)))
	time = dt
	

	for i in range(1,Nt+1):
		u1[1:-1] = -u_1[1:-1] + 2*u[1:-1]+c**2*0.5*\
		((q[1:-1]+q[2:])*(u[2:]-u[1:-1])-(q[1:-1]+q[:-2])*(u[1:-1]-u[:-2])) +\
		dt**2*fval(x[1:-1],t[i],L,w)
		
		if task == "a" or task == "b": #Can be changed to else
			u1[0] = -u_1[0]+2.0*u[0]+ c**2 *(q[0]+q[1])*(u[1]-u[0]) + dt**2*fval(x[0],t[i],L,w)
			u1[-1] = -u_1[-1]+ 2.0*u[-1]+c**2 *(q[-1]+q[-2])*(u[-2]-u[-1]) + dt**2*fval(x[-1],t[i],L,w)

		if task == "c":
			u1[0] = u1[1]
			u1[-1] = u1[-2]

		if task == "d":
			u1[0] = -u_1[0] +2*u[0]+c**2*0.5*(-(q[1]+q[0])*(u[1]-u[0])) + dt**2*fval(x[0],t[i],L,w)
			u1[-1] = -u_1[-1] +2*u[-1]+c**2*0.5*((q[-1]+q[-2])*(u[-1]-u[-2])) + dt**2*fval(x[-1],t[i],L,w)

		u_1[:] = u[:]
		u[:] = u1[:]

	
		if value < np.max(abs(u1[:]-exact(x, (i+1)*dt, L, w)[:])):
			value = np.max(abs(u1[:]-exact(x, (i+1)*dt, L, w)[:]))
			error[:] = u1[:]
			time = (i+1)*dt
		"""
		plt.figure()
		plt.plot(x, exact(x,t[i],L,w),'g')
		plt.plot(x, u1)
		plt.axis([-0.1,2.1,-1.2,1.2])
		plt.xlabel('x')
		plt.ylabel('Values of u')
		plt.title("Exercise %s, numerical solution against analytical \n" \
		 "Time = %d, dt = %.2f, L = %d, dx = %.2f" % (task, T,dt,L,dx))                               	
		plt.savefig('frame_%04d.png' % i)
		"""
	return error, time, x

def convergence_rates(m, w, dx, L, dt, T, V, solver_function, task):
		     
    dt_values = []
    E_values = []
    for i in range(m):
        #u, t = solver_function(V,w,dx,L,dt,T)
        u, t, x = solver(V,w,dx,L,dt,T, task)
        u_e = exact(x, t, L, w)
        E = np.sqrt(dt*np.sum((u_e-u)**2))
        dt_values.append(dt)
        E_values.append(E)
        dt = dt/2

    r = [np.log(E_values[i-1]/E_values[i])/
         np.log(dt_values[i-1]/dt_values[i])
         for i in range(1, m, 1)]
    return r, dt_values, E_values

def main(task):
	import sys
	V = 0.0
	L = 2; T = 6
	dx = 0.09; dt = 0.04
	if dt/dx > 0.5:
		print "dt/dx must be lower than 0.5"
		sys.exit(1)
	w = 1
	m = 5
	#solver(V, w, dx, L, dt, T, task)
	r, Dt, E = convergence_rates(m, w, dx, L, dt, T, V, solver, task) 
	print "Using V = %d, L = %d, dx = %.4f, dt = %.4f" % (V, L, dx, dt)
	print "Dividing dt by 2, %d times" % (m-1)
	for i in range(m-1):
		print "Convergence rate %.5f for dt = %.5f, Error %.5f" % (r[i], Dt[i], E[i])
	

if __name__ == '__main__':
	exercise = ["a", "b", "c", "d"]
	for i in exercise:
		print "Exercise %s" % i 
		main(i)
	
