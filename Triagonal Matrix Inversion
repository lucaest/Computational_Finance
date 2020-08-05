## Triagonal Matrix Inversion

#given data
alpha = np.array([1, 3, 2])
beta = np.array([2, 1])
b = np.array([3, 1, 3])
gamma = np.array([1, 1])
    
def TriagMatrix_Inversion(alpha,beta,gamma,b):
    n = len(alpha)
    alpha_hat = np.zeros(n)
    b_hat = np.zeros(n)
    x = np.zeros(n)

    alpha_hat[0] = alpha[0]
    b_hat[0] = b[0]
    ### bringing the matrix on upper triangular form (forward recursion) according to p.75 in the lecture notes
    for i in range(1,n):
        alpha_hat[i] = alpha[i] - beta[i-1]*gamma[i]/alpha_hat[i-1]
        b_hat[i] = b[i] - b_hat[i-1]*gamma[i]/alpha_hat[i-1]

    ### solving the linear equation system according to (5.17) and (5.18)
    x[n-1] = b_hat[n-1]/alpha_hat[n-1]
    for i in range(n-2,-1,-1):
        x[i] = (b_hat[i]-beta[i] * x[i+1]) / alpha_hat[i]

    return x

#test function
#x = TriagMatrix_Inversion(alpha,beta,gamma,b)

#check function
#x_check = np.linalg.solve([[1, 3, 1], [1, 2, 0], [0, 1, 2]], [1, 3, 3])
