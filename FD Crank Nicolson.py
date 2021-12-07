## Finite Difference Crank Nicolson Scheme for Option Pricing

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

r = 0.05; sigma = 0.2; a = -0.7; b = 0.4; m = 100; nu_max = 2000; T = 1; K = 100

#use matrix inverstion function 

def BS_EuCall_FiDi_CN(r, sigma, a, b, m, nu_max, T, K):
    ### setting the parameters needed for the recursion
    q = 2 * r / sigma ** 2
    dx = (b - a) / m
    dt = sigma ** 2 * T / (2 * nu_max)
    lbd = dt / dx ** 2

    ### range of underlying transformed stock prices
    x = np.arange(a, b + dx, dx)

    ### allocating memory
    w = np.zeros((m + 1, nu_max + 1))

    ### initial values equivalent to transformed payoff at maturity
    w[:, 0] = np.maximum(0, np.exp(x / 2 * (q + 1)) - np.exp(x / 2 * (q - 1)))

    ### setting the main and side diagonals of the tridiagonal matrix 'A_impl'
    alpha = np.ones(m-1) * (1+lbd)
    beta = np.ones(m-1)  * (-lbd/2)
    gamma = np.ones(m-1) * (-lbd/2)

    ### loop over columns of matrix/time
    for i in range(1, nu_max + 1):

        ### loop over rows/underlying (transformed) stock price
        ### note that we do not change the top and bottom row which are equal to zero all the time (simplified boundary conditions)
        for j in range(1, m):
            ### calculating the right hand side of (5.21), can be seen as the explicit part of the CN-scheme
            w[j, i] = lbd/2 * w[j - 1, i - 1] + (1-lbd) * w[j, i - 1] + lbd/2 * w[j + 1, i - 1]

        ### boundary condition for right hand side (next time step)
        #w[m-1,i] = w[m-1,i] + lbd/2 * (np.exp((q + 1) / 2 * b + (q + 1) ** 2 / 4 * i * dt) - np.exp((q - 1) / 2 * b + (q - 1) ** 2 / 4 * i * dt))
        w[1:-1, i] = TriagMatrix_Inversion(alpha, beta, gamma, w[1:-1, i])

        ### boundary condition for explicit part of next iteration (could also be part of line 62) as in lecture notes
        #w[m, i] = np.exp((q + 1) / 2 * b + (q + 1) ** 2 / 4 * i * dt) - np.exp((q - 1) / 2 * b + (q - 1) ** 2 / 4 * i * dt)

    ### retransfoming underlying stock prices
    S = K * np.exp(x[1:-1])

    ### transforming the solution of (5.1) into option prices
    V = K * w[1:-1, nu_max] * np.exp(-x[1:-1] / 2 * (q - 1) - sigma ** 2 / 2 * T * ((q - 1) ** 2 / 4 + q))
    return [S, V]

[S, Vcn] = BS_EuCall_FiDi_CN(r, sigma, a, b, m, nu_max, T, K)
