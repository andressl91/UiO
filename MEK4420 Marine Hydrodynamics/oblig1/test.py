import numpy as np
import matplotlib.pyplot as plt 
n = 5
N = (n-1)*4

#print int(N*0.11),
b = 1.
a = 0.
m = b/2
line = int(N/4.)
print line
x = np.linspace(0,b,line)
x2 = np.zeros(line)
x_vec = np.zeros(N)
y_vec = np.zeros(N)

x2[:] = a + (b-a)*((x[:]-a)/(b-a))**0.2
x2_inv = x2[::-1]
#test = 1-x2[N/2:]
#x2[:N/2] = test[::-1]
x_vec[0] = 0;x_vec[line] = b 
x_vec[2*line] = b ;x_vec[3*line] = 0 
x_vec[:line] = x2[:]; x_vec[line+1:2*line-1] = x2[-1]
x_vec[2*line-1:3*line-2] = x2_inv[:-1]; x_vec[3*line:4*line] = x2[0]

#x_vec[0] = 0;x_vec[line] = b 
#x_vec[2*line] = b ;x_vec[3*line] = 0 
y_vec[:line] = x2[0] ; y_vec[line-1:2*line-1] = x2[:]
y_vec[2*line-1:3*line-2] = x[-1]#; y_vec[3*line:4*line] = x2[0]

print x_vec
plt.figure(1)
plt.plot(x_vec,'*')
#plt.show()



