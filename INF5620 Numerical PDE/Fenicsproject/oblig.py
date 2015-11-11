from dolfin import *
import numpy as np

class solver(object):

    def __init__(self, **kwargs):
    	for key in kwargs:
    		setattr(self, key, kwargs[key])
        self.mesh = UnitSquareMesh(6, 4)
        self.u0 = Expression("1", t=0)
        self.u = None
        self.error = None
	
    def define_mesh(self, mesh):
        self.mesh = mesh

    def define_init(self, exp):
        self.u0 = exp #expression

    def alpha(self):
        return Constant(1)

    def f(self):
        return Constant(0)

    def get_solution(self):
        return self.u

    def get_error(self):
        return self.error
 
    def Piccard(self, iterations):
    	#Define Functionspace
        mesh = self.mesh
    	V = FunctionSpace(mesh, 'CG', 1)
    	u = TrialFunction(V)
    	v = TestFunction(V)

		#Initial Condition    	
    	u0 = self.u0

    	u_1 = interpolate(u0, V)
    	#Define variatonal problem
        dt = self.dt; q = self.q; T = self.T
    	a = u*v*dx + dt/q*self.alpha()*inner(nabla_grad(u), nabla_grad(v))*dx
    	L = (u_1 + dt/q*self.f())*v*dx
    	A = assemble(a)
    	u = Function(V)
    	
        Nt = int(round(T/dt)+1)
        E_list = np.zeros(Nt+1)
    	error = 1
    	k = 0

        ue = interpolate(u0, V)
        e = ue.vector().array() - \
                        ue.vector().array()

        for i in range(1, Nt+1):
            a = 2
            while error > 1E-5 and k < iterations:

                b = assemble(L)
                u0.t = i*dt
                ue = interpolate(u0, V)
                #bc.apply(A, b)
                solve(A, u.vector(), b)
                k += 1
                u_1.assign(u) 
                error = errornorm(u,u_1)

            e = np.max(u_1.vector().array() - \
                        ue.vector().array())
            E = np.sqrt(np.sum(e**2)/ue.vector().array().size)
        self.u = u; self.error = e
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
            #ue = Expression("exp(-pi**2*t)*cos(pi*x", t=0)
            u0 = Expression("cos(pi*x[0])", t=0)
            sol.define_init(u0)
            h = sol.dt
            k = 1
            for i in range(4): 
                dt = sol.dt/2.
                k = 2*k
                Dim2 = UnitSquareMesh(k, k)
                sol.define_mesh(Dim2)
                sol.Piccard(1)
                print sol.error

if __name__ == "__main__":
    #dt = 0.1; q = 0.1; t = 0
    #T = 2.0; test = 10
    sol = solver(dt = 0.1, T = 2.0, q = 0.1)
 
    v = Verify(sol)
    v.analytic()
    #v.analytic()

