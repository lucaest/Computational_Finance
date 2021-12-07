## Simulation of GBM with Exact/Euler/Milstein Method

N = 10000; mu = 0.1; sigma = 0.3; T = 1; X0 = 0.3**2

def Sim_Paths_GeoBM(X0, mu, sigma, T, N):
    
    dt = T/N
    
    Y_exact = X0 * np.ones((N+1))
    Y_euler = X0 * np.ones((N + 1))
    Y_milshtein = X0 * np.ones((N + 1))

    Z = np.random.normal(0, 1, (N+1))

    for i in range(0, N): 
        dW = Z[i] * math.sqrt(dt)
        Y_exact[i] = Y_exact[i-1] * np.exp((mu- math.pow(sigma,2)/2)*dt + sigma * dW)
        Y_euler[i] = Y_euler[i-1] * (1 + mu * dt + sigma * dW)
        Y_milshtein[i] = Y_milshtein[i-1] * (1+mu*dt + sigma*dW + math.pow(sigma,2)/2*(math.pow((dW), 2)- dt))

    return Y_exact, Y_euler, Y_milshtein

samples = Sim_Paths_GeoBM(X0, mu, sigma, T, N)

w = np.arange(0, N+1)*T/N

plt.plot(w, samples[0], linewidth=0.2)
plt.plot(w, samples[1], linewidth=0.2)
plt.plot(w, samples[2], linewidth=0.2)

