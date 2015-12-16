import numpy as np

a = np.zeros(4)
a[1] = -1.2
a[0] = 1

i = np.where(a<-1)
a[np.where(a<-1)] = -1
print a