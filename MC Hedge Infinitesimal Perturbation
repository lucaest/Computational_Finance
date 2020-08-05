## MC Hedge with Infinitesimal Perturbation

t = 0; S0 = 100; r = 0.05; sigma = 0.2; T = 1; N = 10000; K = 90

def vanilla_call_payoff(x):
    return max(x - K, 0)

def EuOptionHedge_BS_MC_IP (St, r, sigma, g, T, t, N):
    
    def derivative(x):
        return scipy.misc.derivative(g,St * np.exp((r-sigma**2/2)*(T-t)+sigma*np.sqrt(T-t)*x))
    
    delta = np.zeros((N))
    X = np.random.normal(0,1,N)
    
    for i in range(0,N):
        delta[i] = np.exp(-(np.power(sigma,2)/2)*(T-t)+sigma*np.sqrt(T-t)*X[i]) * derivative(X[i])

    return np.mean(delta)

phi1 = EuOptionHedge_BS_MC_IP(S0, r, sigma, vanilla_call_payoff, T, t, N)

#Hedge according to p.48 in lecture notes
def EuCallHedge_BlackScholes(t, S_t, r, sigma, T, K):
    d_1 = (math.log(S_t / K) + (r + 1 / 2 * math.pow(sigma, 2)) * (T - t)) / (sigma * math.sqrt(T - t))
    return scipy.stats.norm.cdf(d_1)

phi1_bs = EuCallHedge_BlackScholes(t, S0, r, sigma, T, K)
