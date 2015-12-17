import numpy as np

class Figure:
	def __init__(self, r_a, N):
		self.r_a = r_a; self.N = N; 
		self.x = None; self.y = None
		self.mat = None
	def figure_points(self):
		N = self.N; r_a = self.r_a
		N = N/4 *4
		Np = N/4
		x = np.zeros(N+1)
		y = np.zeros(N+1)
		d1 = -r_a
		d2 = r_a
		for i in range(Np+1):
			x[i] = d1 +(d2-d1)/2 *(1-np.cos(i*np.pi/Np))
			y[i] = -r_a
			x[i+Np] = r_a
			y[i+Np] = d1 +(d2-d1)/2 *(1-np.cos(i*np.pi/Np))
			x[i+2*Np] = -(d1 +(d2-d1)/2 *(1-np.cos(i*np.pi/Np)))
			y[i+2*Np] = r_a
			x[i+3*Np] = -r_a
			y[i+3*Np] = -(d1 +(d2-d1)/2 *(1-np.cos(i*np.pi/Np)))
		self.x = x; self.y = y
		return x ,y

	def angles(self):
		N = self.N; 
		x, y = self.figure_points()
		mat = np.zeros((N,N))
		theta = np.zeros(N)
		for i in range(N):
			initx = (x[i+1] + x[i])*0.5
			inity = (y[i+1] + y[i])*0.5

			test = ((x[:-1]-initx)*(x[1:]-initx) + (y[:-1]-inity)*(y[1:]-inity))\
			/(np.sqrt((x[:-1]-initx)**2+(y[:-1]-inity)**2)\
			*np.sqrt((x[1:]-initx)**2+(y[1:]-inity)**2)) 
			test[np.where(test<-1)] = -1 #Secure valid values
			test[np.where(test>1)] = 1
			theta = -np.arccos(test[:])
			
			mat[i][:i+1] = theta[:(i+1)]
			mat[i][i] = -np.pi
			mat[i][i+1:] = theta[i+1:]
		self.mat = mat
		return mat

	def integrate(self, direction):
		area = 0
		N = self.N
		#x = self.x; y = self.y
		a = self.r_a
		x, y = self.figure_points()
		integral = np.zeros(N)
		x_0 = 0.5*(x[:-1]+x[1:])
		y_0 = 0.5*(y[:-1]+y[1:])

		#x_m = (x[1:]-x[:-1])
		#y_m	= (y[1:]-y[:-1])
		ds = np.sqrt((x[1:]-x[:-1])**2 + (y[1:]-y[:-1])**2)
		if direction == 11:
			#HSUK BYTTE
			n = -(x[:-1]-x[1:])/ds
		if direction == 22:
			n = -(y[:-1]-y[1:])/ds
		if direction == 66:
			ny = (x[1:]-x[:-1])/ds
			nx = -(y[1:]-y[:-1])/ds
			rx = x_0
			ry = y_0
			n = rx*ny - nx*ry	

		for i in range(N): #N-1
			radius_a = np.sqrt((x_0[i] -x[:-1])**2 + (y_0[i] - y[:-1])**2)
			radius_b = np.sqrt((x_0[i] -x[1:])**2 + (y_0[i] - y[1:])**2)
			f_a = np.log(radius_a[:]) 
			f_b = np.log(radius_b[:])
			integral[i] = np.sum(0.5*(f_a[:] + f_b[:])*n[:]*ds[:])
		return integral

	def added(self, sol, direction):
		x = self.x; y = self.y
		a = self.r_a
		x_0 = 0.5*(x[:-1]+x[1:])
		y_0 = 0.5*(y[:-1]+y[1:])
		ds = np.sqrt((x[1:]-x[:-1])**2 + (y[1:]-y[:-1])**2)
		if direction == 11:
			#Pass paa fort
			num = x
			n = -(num[:-1]-num[1:])/ds
		if direction == 22:
			num = y
			n = -(num[:-1]-num[1:])/ds
		if direction == 66:
			ny = (x[1:]-x[:-1])/ds
			nx = -(y[1:]-y[:-1])/ds
			rx = x_0
			ry = y_0
			n = rx*ny - nx*ry	

		area = 0.5*np.sum((sol+sol)*n*ds)
		#area = 0.5*np.sum((sol[:-1]+sol[1:])*n*ds)	
		return area

	def solve(self, direction):
		fig = Figure(self.r_a, self.N)
		mat = fig.angles()
		integ = fig.integrate(direction)
		test = np.linalg.solve(mat, integ)
		return fig.added(test, direction)


if __name__ == "__main__":
	a = 4; N = 100
	fig = Figure(a, N)
	fig.angles()
	#fig.solve(11)
	case = "1"
	if case == "1":
		a = 4; N = 800
		fig = Figure(a, N)
		direction = [11, 22, 66]
		exact = [1.51*np.pi*(a)**2,1.51*np.pi*(a)**2 , 0.234*np.pi*a**4]
		print "-----------------------------------"
		print "Square"
		print "Lenght chosen as %d, with %d nodes" % (a,N)
		print "Added mass is Calculated as"
		for i in range(len(direction)):
			print " For direction %.d Numerical solution %.3f, exact solution %.3f" \
				% (direction[i], fig.solve(direction[i]), exact[i])
			print "Error %.3f %%" % ((np.abs(fig.solve(direction[i])))/exact[i])
		print "-----------------------------------"