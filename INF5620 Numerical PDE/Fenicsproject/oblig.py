from dolfin import *
import numpy as np

class solver(object):

    def __init__(self, **kwargs):
    	for key in kwargs:
    		setattr(self, key, kwargs[key])
        self.mesh = UnitSquareMesh(6, 4)
        self.u0 = None
        self.u = None
        self.f = None
        self.alpha = None
	
    def define_mesh(self, mesh):
        self.mesh = mesh

    def define_init(self, exp):
        self.u0 = exp #expression

    def define_dt(self, dt):
        self.dt = dt

    def define_f(self, f):
        self.f = f

    def define_alpha(self, alpha):
        self.alpha = alpha

    def alpha_func(self, u=0):
        #return 1+ u*u
        return self.alpha(u)

    def f_func(self,t=0):
        return self.f

    def get_solution(self):
        return self.u

    def get_error(self):
        return self.error
 
    def Piccard(self, iterations, getstep = None):
    	#Define Functionspace
        mesh = self.mesh

    	V = FunctionSpace(mesh, 'CG', 1); self.V = V
    	u = TrialFunction(V)
    	v = TestFunction(V)

		#Initial Condition    	
    	u0 = self.u0
    	u_1 = interpolate(u0, V)
    	#Define variatonal problem
        dt = self.dt; q = self.q; T = self.T
    	a = u*v*dx + dt/q*self.alpha_func(u_1)*inner(nabla_grad(u), nabla_grad(v))*dx
    	L = (u_1 + dt/q*self.f_func() )*v*dx
    	A = assemble(a)
    	u = Function(V)
        Nt = int(round(T/dt)+1)
        E_list = np.zeros(Nt+1)
    	error = 1
    	k = 0
 
        for i in range(1, Nt+1):
            #a = 2
            while error > 1E-5 and k < iterations:
                b = assemble(L)
                self.f_func().t = i*dt 
                #bc.apply(A, b)
                solve(A, u.vector(), b)
                k += 1
                u_1.assign(u) 
                error = errornorm(u,u_1)

            if getstep == i*dt:
                u0.t = i*dt
                ue = interpolate(u0, V)
                #print ue.vector().array()
                e = ue.vector().array() - \
                    u_1.vector().array()
                E = np.sqrt(np.sum(e*e)/u_1.vector().array().size)
                self.u = u_1
                return E
        self.u = u #This is only the last
    	return u
		
class Verify(object):

        def __init__(self, solverclass):
            self.sol = solverclass

        def constant(self):
            dim = 1
            Dim1 = IntervalMesh(10,0,6)
            Dim2 = UnitSquareMesh(6, 4)
            Dim3 = UnitCubeMesh(6, 10, 5)
            Dim = [Dim1, Dim2, Dim3]
            #sol.define_alpha(Constant(1))
            sol.define_alpha(lambda u: Constant(1))
            sol.define_f(Constant(0))
            sol.define_init(Expression("1", t=0))
            while dim <= 3:
                sol.define_mesh(Dim[dim-1])
                u = self.sol.Piccard(1)
                u_vector = u.vector().array()
                #Have to use around due to small errors in solver
                u_vector = np.around(u_vector)
                test = all(x == u_vector[0] for x in u_vector)
                assert(test == True)
                dim += 1

        #Checking convergencerate for 
        #u_e = exp(-pi**2*t)cos(pi*x)
        def analytic(self):
            ue = Expression("exp(-pi*pi*t)*cos(pi*x[0])", t=0)
            #u0 = Expression("cos(pi*x[0])", t=0)
            sol.define_init(ue)
            #sol.define_alpha(Constant(1))
            sol.define_alpha(lambda u: Constant(1))
            sol.define_f(Constant(0))
            dt = sol.dt
            Nx = int(1./np.sqrt(dt))
            getstep = dt*3
            for i in range(5): 
                Dim2 = UnitSquareMesh(Nx, Nx)
                sol.define_mesh(Dim2)
                error = sol.Piccard(1, getstep)
                print error/dt
                dt = sol.dt/2.
                Nx = int(1./np.sqrt(dt))
                sol.define_dt(dt)


        def manufactured(self):
            Dim1 = IntervalMesh(1000,0,1)
            sol.define_mesh(Dim1)
            ue = Expression("t*x[0]*x[0]*(0.5-x[0]/3.0)", t=0)
            sol.define_init(ue)
            sol.define_f(Expression("-q*x[0]*x[0]*x[0]/3\
                            + q*x[0]*x[0]/2 + \
            8*t*t*t *x[0]*x[0]*x[0]*x[0]*x[0]*x[0]*x[0]/9 -\
            28*t*t*t *x[0]*x[0]*x[0]*x[0]*x[0]*x[0]/9 + \
            7*t*t*t *x[0]*x[0]*x[0]*x[0]*x[0]/2 - \
            5*t*t*t *x[0]*x[0]*x[0]*x[0]/4 + \
            2*t*x[0]-t", t=0,q=sol.q))
            sol.define_alpha(lambda u: 1+u*u)
            #getstep = [i*sol.dt for i in [1,3]]
            getstep = [.1]
            for i in getstep:
                sol.Piccard(1, i)
                ue.t = i
                ue = interpolate(ue, sol.V)
                plot(sol.u)
                plot(ue)
                interactive()
                #import matplotlib.pyplot as plt
                #plt.figure(1)
                #plt.plot(sol.u.vector().array())
                #plt.plot(ue.vector().array(), 'r')
                #plt.show()
                

if __name__ == "__main__":
    #dt = 0.1; q = 0.1; t = 0
    #T = 2.0; test = 10
    sol = solver(dt = 0.01, T = 1.0, q = 1.0)
 
    v = Verify(sol)
    #v.constant()
    #v.analytic()
    v.manufactured()
