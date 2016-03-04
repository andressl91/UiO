#Defining boundary for neumann condition
mf = MeshFunction("size_t", mesh, 1) #3rd argument is dim of egde
mf.set_all(0)

nmann = Neumann()
nmann.mark(mf, 1)
#Define new ds based of subdomain MeshFunction
ds = Measure("ds")[mf]

class Neumann(SubDomain):
    def inside(self, x, on_boundary):
        return on_boundary and ( near(x[1], 1) or near(x[1], 0) )
