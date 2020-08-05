## MC European Call

S0 = 100; r = 0.05; sigma = 0.2; T = 1; M = 10000; K = 90

def eu_call_payoff(x):
    return max(x - K, 0)

def Eu_Option_BS_MC (S0, r, sigma, T, K, M, f):
    
    S = np.zeros((M+1))
    V = np.zeros((M+1))
    X = np.random.normal(0, 1, M+1)
    
    for i in range(0, M+1):
        S[i] = S0 * math.exp((r-sigma**2/2)*T + sigma*math.sqrt(T)*X[i])
        V[i] = eu_call_payoff(S[i])    
    
    V0 = math.exp(-r*T) * np.mean(V)
    
    epsilon = 1.96 * math.sqrt(np.var(V) / M)
    c1 = V[0] - epsilon
    c2 = V[0] + epsilon
    
    return V0, c1, c2

V = Eu_Option_BS_MC(S0, r, sigma, T, K, M, eu_call_payoff)

def bs_formula(S0, T, K, r, sigma):
    
    d1 = (math.log(S0/K) + r*T + 0.5*sigma**2*T)/sigma*T
    d2 = (math.log(S0/K) + r*T - 0.5*sigma**2*T)/sigma*T
    V = S0 * scipy.stats.norm.cdf(d1) - K*math.exp(-r*T)*scipy.stats.norm.cdf(d2)
   
    return V

V_bs = bs_formula(S0, T, K, r, sigma)
