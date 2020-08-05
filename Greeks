## Greeks approximation

S0 = 120
#S0 = np.arange(60, 140+1); 
r = 0.03; sigma = 0.2; K = 100; T = 1; eps = 0.001

def vanilla_call_payoff(x):
    return max(x - K, 0)

# compute Black-Scholes price by integration
def BS_Price_Int(S0, r, sigma, T, f):
    # define integrand as given in the exercise
    def integrand(x):
        return 1 / math.sqrt(2 * math.pi) * f(
            S0 * math.exp((r - 0.5 * math.pow(sigma, 2)) * T + sigma * math.sqrt(T) * x)) * math.exp(-r * T) * math.exp(
            -1 / 2 * math.pow(x, 2))

    # perform integration
    I = integrate.quad(integrand, -np.inf, np.inf)
    # return value of the integration
    return I[0]

def BS_Greeks_num(r, sigma, S0, T, g ,eps):
    
    V0 = BS_Price_Int(S0, r, sigma, T, g)
    
    Vplus_delta = BS_Price_Int(S0+eps*S0, r, sigma, T, g)
    delta = (Vplus_delta - V0)/(eps*S0)
    
    Vplus_vega = BS_Price_Int(S0, r, sigma+eps*sigma, T, g)
    vega = (Vplus_vega - V0)/(eps*sigma)
    
    Vplus_gamma = BS_Price_Int(S0+eps*S0, r, sigma, T, g)
    Vminus_gamma = BS_Price_Int(S0-eps*S0, r, sigma, T, g)
    gamma = (Vplus_gamma - 2*V0 + Vminus_gamma)/((eps*S0)**2)
    
    return delta, vega, gamma

#[delta, vega, gamma] = BS_Greeks_num(r, sigma, S0, T, vanilla_call_payoff, eps)

delta = np.zeros((N))
gamma = np.zeros((N))
vega = np.zeros((N))

for i in range(0, len(S0)):
    delta[i] = BS_Greeks_num(r, sigma, S0[i], T, vanilla_call_payoff, eps)[0]
    vega[i] = BS_Greeks_num(r, sigma, S0[i], T, vanilla_call_payoff, eps)[1]
    gamma[i] = BS_Greeks_num(r, sigma, S0[i], T, vanilla_call_payoff, eps)[2]

plt.plot(S0, delta)
plt.plot(S0, vega)
plt.plot(S0, gamma)
