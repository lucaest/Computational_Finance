## CRR Option Pricing

S0 = 100; r = 0.05; sigma = 0.3; T = 1; #M = np.arange(10, 500+1); 
M = 100
K = 120; 

def eu_put_payoff(x):
    return np.maximum(K - x, 0)

def CRR_AmEuPut (S_0, r, sigma, T, M, K, eu):
    
    dt = T/M
    beta = 0.5 * (math.exp(-r*dt) + math.exp((r+sigma**2)*dt))
    u = beta + math.sqrt(beta**2 - 1)
    d = 1/u
    q = (math.exp(r*dt) - d) / (u-d)
    
    S = np.zeros((M+1, M+1))
    V = np.zeros((M+1, M+1))
    
    S[0, 0] = S0
    
    
    for i in range(0, M+1):
        for j in range(0, i+1):
            S[j, i] = S[0, 0] * u**j * d**(i-j)
    
    V[:, M] = eu_put_payoff(S[:, M])

    for i in range(M-1, -1, -1):
        for j in range(0, i+1):
            if eu == 1:#european put
                V[j, i] = np.exp(-r*dt) * (q*V[j+1, i+1] + (1- q)*V[j, i+1])
            if eu == 0:#american put
                V[j, i] = np.maximum(eu_put_payoff(S[j, i]), np.exp(-r*dt)*(q*V[j+1, i+1] + (1- q)*V[j, i+1]))
           
    return V[0, 0]

print(CRR_AmEuPut(S0, r, sigma, T, M, K, 0))

