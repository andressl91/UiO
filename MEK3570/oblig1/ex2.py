import numpy as np

#Dimention of matrix
N = 50

#Construct random Matrix of dimention N
A = np.random.randint(5, size=(N, N))

################### 1.94 ###################

#Construct identity Matrix of dimention N
I = np.eye(N)

#Calculation of I : A 
I_leftside = 0
for i in range(N):
    for j in range(N):
        I_leftside += I[i, j]*A[i, j]

#Calculation of trace of matrix A
A_trace = float(np.trace(A))

#Calculation of A : I
I_rightside = 0
for i in range(N):
    for j in range(N):
        I_rightside += A[i, j]*I[i, j]

test_1 = I_leftside == A_trace
test_2 = A_trace == I_rightside
print("Property 1.94")
print("I : A = tr(A) is %r" % test_1)
print("tr(A) = A : I is %r \n" % test_2)  

################### 1.96 ###################

u = np.random.randint(100, size=N)
v = np.random.randint(100, size=N)

u_cross_v = np.outer(u, v)
A_leftside = 0
A_rightside = 0
for i in range(N):
    for j in range(N):
        A_leftside  += A[i, j]*u_cross_v[i, j]
        A_rightside += u_cross_v[i, j]*A[i, j]

u_Av = np.dot(u, np.dot(A,v))

test_1 = A_leftside == u_Av
test_2 = u_Av == A_rightside

print("Property 1.96")
print("A : (u (x) v) = u dot Av is %r" % test_1)
print("u dot Av = (u (x) v) : A is %r" % test_2)  




