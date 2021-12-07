## MC Option Pricing with Importance Sampling


S0 = 100; r = 0.05; sigma = 0.3; K = 220; T = 1; N = 10000
# compute optimal mu
d = ((np.log(K/S0)) - (r - pow(sigma, 2)/2) * T) / sigma * math.sqrt(T)
mu = d

def BS_EuCall_MC_IS (S0, r, sigma, K, T, mu, N):
    
    def eu_call_payoff(x):
        return max(x - K, 0)
    
    S1 = np.zeros((N+1))
    S2 = np.zeros((N+1))
    V = np.zeros((N+1))
    V_IS = np.zeros((N+1))
    S1[0] = S0
    S2[0] = S0
    X = np.random.normal(0, 1, N+1)
    Y = np.random.normal(mu, 1, N+1)
    
    for i in range(1, N+1):
        S1[i] = S0 * math.exp((r-sigma**2/2)*T + sigma*math.sqrt(T)*X[i])
        S2[i] = S0 * math.exp((r-sigma**2/2)*T + sigma*math.sqrt(T)*Y[i])
        V[i] = math.exp(-r*T) * eu_call_payoff(S1[i])
        V_IS[i] = math.exp(-r * T-Y[i]*mu+(mu**2)/2) * eu_call_payoff(S2[i])
    
    V0 = np.mean(V)
    V_IS0 = np.mean(V_IS)
        
    return V0, V_IS0

np.random.seed(123)
Vs = BS_EuCall_MC_IS(S0, r, sigma, K, T, mu, N)
