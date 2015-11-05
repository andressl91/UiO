# http://hplgit.github.io/decay-book/doc/pub/book/sphinx/._book007.html#decay-app-nuclear
# Given lambda_, dt, N

import numpy as np

dt = 0.1
lambda_ = 1
N = 100
T = 10.0
dt = T/N + 1


def decay(dt, N, lambda_):
	#Decay takes the parameters T for time, N for number of nucelus
	#lambda_ as parameter for probability. Calculates the prob for
	#decay and calc number of N who decayed
	#dt = T/(N+1)

	p = lambda_*dt #Probability for decay
	outcome = np.random.binomial(N, p)
	return N-outcome


def main():
	count = 0
	lambda_ = 1.0
	N = 100
	T = 10.0
	dt = T/(N+1)
	while dt <= T:
		res = decay(dt, N, lambda_)
		N -= res
		dt = T/(N+1)
		count += 1
		

	

if __name__ == '__main__':
	main()


