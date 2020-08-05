## MC Option Pricing with Control Variate


S0 = 100; r = 0.05; sigma = 0.3; K = 110; T = 1; M = 100000

def BS_EuOption_MC_CV (S0, r, sigma, T, K, M):
    
    def quanto_call_payoff(x):
        return max((x - K), 0)*x
    
    def vanilla_call_payoff(x):
        return max(x - K, 0)
    
    def bs_formula(S0, T, K, r, sigma):
        
        d1 = (math.log(S0/K) + r*T + 0.5*sigma**2*T)/sigma*T
        d2 = (math.log(S0/K) + r*T - 0.5*sigma**2*T)/sigma*T
        V = S0 * scipy.stats.norm.cdf(d1) - K*math.exp(-r*T)*scipy.stats.norm.cdf(d2)
       
        return V

    S1 = np.zeros((M+1))
    #S2 = np.zeros((M+1))
    #V1 = np.zeros((M+1))
    #V2 = np.zeros((M+1))
    V_CV = np.zeros((M+1))
    X = np.random.normal(0, 1, M+1)
    #Y = np.random.normal(0, 1, M+1)  
    
    S1[0] = S0
    #S2[0] = S0
    
    beta = 0.5
    
    for i in range(1, M+1):
        S1[i] = S0 * math.exp((r-sigma**2/2)*T + sigma*math.sqrt(T)*X[i])
        #S2[i] = S0 * math.exp((r-sigma**2/2)*T + sigma*math.sqrt(T)*Y[i])   
        #V1[i] = quanto_call_payoff(S1[i])
        #V2[i] = vanilla_call_payoff(S1[i])
        V_CV[i] = quanto_call_payoff(S1[i])- beta * vanilla_call_payoff(S1[i])
    
    V_CV_0 = math.exp(-r*T) * np.mean(V_CV) + beta * bs_formula(S0, T, K, r, sigma)
    
    return V_CV_0

np.random.seed(123)
price = BS_EuOption_MC_CV(S0, r, sigma, T, K, M)
