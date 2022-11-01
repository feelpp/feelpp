import numpy as np
#from scipy import linalg


def first_difference(T, Y):
	"""
	T is a vector of two times, greatest to least
	Y is a vector older differences, [dj,d(j+1)]
	See algorthm in paper
	"""
	return (Y[0]-Y[1])/(T[0] - T[1])


def backward_differences(T):
	"""
	Generate the divided difference coefficients. T is a vector of 
	times from least to greatest
	T = [t_n,t_n+1,....,t_m+m]
	"""
	numOfTimes = len(T)
	#the number of steps in the method
	m = numOfTimes - 1
	#generate the initial differences, which
	#is just the standard basis.
	D = np.array([[np.float64((i+1) == (numOfTimes-j))
	             for i in range(numOfTimes)] for j in range(numOfTimes)])
	differences = np.zeros_like(D)
	differences[0] = D[0]

	for q in range(1, numOfTimes):
		for j in range(numOfTimes - q):
			D[j] = first_difference([T[m-j], T[m-j-q]], [D[j], D[j+1]])
			differences[q] = D[0]
	return differences


def bdf_coefficients(T):
	differences = backward_differences(T)
	m = len(T)-1
	return np.sum(np.prod([T[m]-T[m-i] for i in range(1, j)])*differences[j] for j in range(1, m+1))


def bdf_coefficients_and_differences(T, order):
	differences = backward_differences(T)
	m = len(T)-1
	#calculate filter coefficient for increasing order by one
	eta = np.prod([T[m]-T[m-i] for i in range(1, m)]) / \
            np.sum(1./(T[m] - T[m-j]) for j in range(1, m+1))

	return [np.sum(np.prod([T[m]-T[m-i] for i in range(1, j)])*differences[j] for j in range(1, order+1)), differences, eta]


print(backward_differences([1,2,3,4,5]))

T = [1,2,3]
print(backward_differences(T))

m = 4
#for j in range(1,m+2):
#	print(1/(np.prod([T[m]-T[m-i] for i in range(1,j)])))

#differences([1.0,2.0,3.0,4.5])


#backward_differences([1.0,2.0,3.0,4.0,5.5])
[c,d,eta]=bdf_coefficients_and_differences([1, 2, 3,4], 3)
print("== bdf c: ", c)
print("== bdf d: ", d)
print("== bdf eta: ", eta)

#print(bdf_coefficients([1,2,3.1]))
