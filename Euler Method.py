## Simulation of PDE paths with Euler Method

X0 = 0.3**2
gamma0 = 0.3**2; sigma = 0.2; lbd = 2.5; kappa = 0.3**2; N = 10000; m = 100; T = 1

def a(kappa, lbd, gamma):
    return kappa - lbd*gamma

def b(gamma, sigma):
    return np.sqrt(gamma) * sigma
        

def SimPath_Ito_Euler(X0, a, b, T, m, N):
    
    dt = T/m
    Y = np.zeros((N+1, m+1))
    Y[:, 0] = X0
    Z = np.random.normal(0, 1, (N+1, m+1))
    
    for i in range(0, N):    
        for j in range(1, m+1):
            dW = Z[i, j] * math.sqrt(dt)
            Y[i, j] = Y[i, j-1] + a(kappa, lbd, Y[i, j-1])*dt + b(Y[i, j-1], sigma)*dW

    return Y

samples = SimPath_Ito_Euler(gamma0, a, b, T, m, N)

#plt.plot(samples[75, :] , linewidth=0.5)
for i in range(75, 80):
    plt.plot(samples[i, :] , linewidth=0.5)
