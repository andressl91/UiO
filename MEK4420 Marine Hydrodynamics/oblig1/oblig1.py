import numpy as np
import matplotlib.pyplot as plt

class Figure:
	def __init__(self, r_a, r_b, N):
		self.r_a = r_a ;self.r_b = r_b
		self.N = N; self.phi = None
		self.x = None; self.y = None
		self.mat = None
	def figure_points(self):
		N = self.N; r_a = self.r_a ; r_b = self.r_b 
		dx = 2*np.pi/(N)
		x = np.zeros(N+1)
		y = np.zeros(N+1)
		for i in range(N+1):
			x[i] = r_a*np.cos(dx*i)
			y[i] = r_b*np.sin(dx*i)
			self.x = x; self.y = y
		return x,y

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

			theta = -np.arccos(test[:])
			
			mat[i][:i+1] = theta[:(i+1)]
			mat[i][i] = -np.pi
			mat[i][i+1:] = theta[i+1:]
		self.mat = mat
		return mat
	def integrate(self, direction):
		area = 0
		N = self.N
		a = self.r_a; b = self.r_b
		x, y = self.figure_points()
		integral = np.zeros(N)
		x_0 = 0.5*(x[:-1]+x[1:])
		y_0 = 0.5*(y[:-1]+y[1:])
		ds = np.sqrt((x[1:]-x[:-1])**2 + (y[1:]-y[:-1])**2)
		if direction == 11:
			num = x 
			n = -(num[:-1]+num[1:])/(2.*a**2)/np.sqrt(((x[:-1]+x[1:])/(2.*a**2))**2 + ((y[:-1]+y[1:])/(2.*b**2))**2)
		if direction == 22:
			num = y
			n = -(num[:-1]+num[1:])/(2.*b**2)/np.sqrt(((x[:-1]+x[1:])/(2.*a**2))**2 + ((y[:-1]+y[1:])/(2.*b**2))**2)
		if direction == 66:
			nx = (x_0/(a**2))/np.sqrt((x_0**2/(a**4)) + ((y_0**2)/(b**4)) )
			ny = (y_0/(b**2))/np.sqrt((x_0**2/(a**4)) + ((y_0**2)/(b**4)) )
			rx = x_0
			ry = y_0
			n = x_0*ny - y_0*nx
		
		for i in range(N): #N-1
			radius_a = np.sqrt((x_0[i] -x[:-1])**2 + (y_0[i] - y[:-1])**2)
			radius_b = np.sqrt((x_0[i] -x[1:])**2 + (y_0[i] - y[1:])**2)
			f_a = np.log(radius_a[:]) 
			f_b = np.log(radius_b[:])
			integral[i] = np.sum(0.5*(f_a[:] + f_b[:])*n[:]*ds[:])
		return integral

	def added(self, sol, direction):
		x = self.x; y = self.y
		a = self.r_a; b = self.r_b
		area = 0
		ds = np.sqrt((x[1:]-x[:-1])**2 + (y[1:]-y[:-1])**2)
		if direction == 11:
			num = x
			n = -(num[:-1]+num[1:])/(2.*a**2)/np.sqrt(((x[:-1]+x[1:])/(2.*a**2))**2 + ((y[:-1]+y[1:])/(2.*b**2))**2)
		if direction == 22:
			num = y
			n = -(num[:-1]+num[1:])/(2.*b**2)/np.sqrt(((x[:-1]+x[1:])/(2.*a**2))**2 + ((y[:-1]+y[1:])/(2.*b**2))**2)
		if direction == 66:
			x_0 = 0.5*(x[:-1]+x[1:])
			y_0 = 0.5*(y[:-1]+y[1:])
			nx = -(x_0/(b**2)) / np.sqrt((x_0**2/(a**4)) + ((y_0**2)/(b**4)) )
			ny = -(y_0/(a**2)) / np.sqrt((x_0**2/(a**4)) + ((y_0**2)/(b**4)))
			rx = x_0
			ry = y_0
			n = ny*rx - nx*ry

		area = np.sum((sol)*n*ds)	
		return area

	def solve(self, direction):
		fig = Figure(self.r_a,self.r_b, self.N)
		mat = fig.angles()
		integ = fig.integrate(direction)
		test = np.linalg.solve(mat, integ)
		self.phi = test
		return fig.added(test, direction)
		
if __name__ == "__main__":
	for j in [1,2]:
		r_a = 2; r_b = 2; N = 200*j
		angle = [i*2.*np.pi/N for i in range(N)]
		plt.figure(1)

		fig = Figure(r_a,r_b, N)
		direction = [11, 22, 66]
		exact = [np.pi*r_a**2, np.pi*r_b**2, 10**-10]
		print "---------------------------BEGIN SIMULATION-----------------------------------------------------"
		print
		print "--------------------------- -----CIRCLE-----------------------------------------------------"
		print "Radius chosen as %d, with %d nodes" % (r_a,N)
		for i in range(len(exact)):
			print "--------------------------------------------------------------------------------"
			print " For direction %.d Numerical solution %.3f, exact solution %.3f, error %.3f %%" \
				% (direction[i], fig.solve(direction[i]), exact[i],(abs(fig.solve(direction[i])))/exact[i])
			if i == 0 and N== 200:
				larg_error = np.max(fig.phi-(-r_a*np.cos(angle)))
				plt.plot(fig.phi, label = 'Numerical')
				plt.plot(-r_a*np.cos(angle),label = 'Exact')
				plt.xlabel("Node")
				plt.title("Potential on the Circle for radius = %d, N = %d, \n \
Calculated added mass %.4f, Exact solution added mass %.4f" % (r_a, N, fig.solve(direction[0]),exact[0]))
				plt.legend(loc='upper left')
				plt.savefig('present1.png')
		print

		print "--------------------------------Ellipse-----------------------------------------------------"
		r_a = 1; r_b = 3; 
		fig = Figure(r_a,r_b, N)
		direction = [11, 22, 66]
		exact = [np.pi*r_b**2, np.pi*r_a**2, 1/8.*np.pi*(r_a**2-r_b**2)**2]
		print "Radius r_a = %d, r_b = %d , with %d nodes" % (r_a,r_b,N)
		for i in range(len(exact)):
			print "--------------------------------------------------------------------------------"
			print " For direction %.d Numerical solution %.3f, exact solution %.3f" \
				% (direction[i], fig.solve(direction[i]), exact[i])
			print "Error %.3f %%" % ((abs(fig.solve(direction[i])-exact[i]))/exact[i])

		