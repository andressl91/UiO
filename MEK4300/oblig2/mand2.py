from dolfin import *

#mesh = Mesh(tube.xml)
mesh = UnitSquareMesh(4,4) 
print "Tittei"
plot(mesh)
interactive()
