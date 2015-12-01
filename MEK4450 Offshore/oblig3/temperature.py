import numpy as np 

class Pipe:

	def __init__(self, **kwargs):
		self.wellhead = []
		self.shore = [] 
		self.wellpipe = []
		self.shorepipe = []
		for key in kwargs:
			setattr(self, key, kwargs[key])
			

	def temp(self, x, Q):
		T = self.Tsea + (self.Tres - self.Tsea) \
			*np.exp(-(np.pi*self.U*self.D*x/ \
					(Q*self.Cp)))
		return T

	def pipe2(self):
		self.D = 12*0.0254 #Pipe diameter in meters
		self.Tsea = 7. #Average Temperature sea and reservoir
		self.U = 20. #Heat transfer

	def plot(self, Q, x, T, show = True):

		import matplotlib.pyplot as plt
		for i in range(len(Q)):
			fig1 = plt.figure(1)
			ax1 = fig1.add_subplot(111)
			ax1.plot(x, self.temp(x,Q[i]))
			self.wellhead.append(self.temp(x[-1],Q[i]))
			self.wellpipe.append(self.temp(x,Q[i]))

		ax1.set_xlabel("Distance meters")
		ax1.set_ylabel("Temperature Celcius")	
		ax1.legend(["Q = %dKg/s" % b for b in Q ])
		#plt.savefig("wellhead.png")

		welltemp = [pipe.temp(x[-1], Q[:])]
		#From wellhead to shore

		self.pipe2()
		welltemp = np.asarray(welltemp)
		Q = np.asarray([60, 50, 40, 34, 28])
		x = np.linspace(0, 50000, 50001)
		for j in range(5):
			self.Tres = welltemp[0,j]
			fig2 = plt.figure(2)
			ax2 = fig2.add_subplot(111)
			ax2.plot(x, self.temp(x,Q[j]))
			self.shore.append(self.temp(x[-1],Q[i]))
			self.shorepipe.append(self.temp(x,Q[i]))

		ax2.set_xlabel("Distance meters")
		ax2.set_ylabel("Temperature Celcius")	
		ax2.legend(["Q = %dKg/s" % b for b in Q ])
		#plt.savefig("shore.png")
		if show == True:

			plt.show()

#From reservoir to well

Q = np.asarray([30, 25, 20, 17, 14])
x = np.linspace(0,4000,4001)
#From reservoir to wellhead
pipe = Pipe(
		Tres = 150., #Temperature reservoir
		D = 5*0.0254, #Pipe diameter in meters
		Tsea = (150.+7.)/2., #Average Temperature sea and reservoir
		U = 9., #Heat transfe
		Cp = 2000., #Heat capasity
		)
welltemp = [pipe.temp(x[-1], Q[:])]
pipe.plot(Q,x, welltemp)



