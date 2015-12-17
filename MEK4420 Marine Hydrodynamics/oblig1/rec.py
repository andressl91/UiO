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

			theta = -np.arccos(test[:])
			
			mat[i][:i+1] = theta[:(i+1)]
			mat[i][i] = -np.pi
			mat[i][i+1:] = theta[i+1:]
		self.mat = mat
		return mat

a = 1; N = 8
fig = Figure(a, N)