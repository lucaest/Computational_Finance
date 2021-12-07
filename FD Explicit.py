## Finite Difference Explicit Scheme for Option pricing

r = 0.05; sigma = 0.2; a = -0.7; b = 0.4; m = 100; nu_max = 2000; T = 1; K = 100

def BS_EuCall_FiDi_Explicit(r, sigma, a, b, m, nu_max, T, K):
    ### setting the parameters needed for the recursion
    q = 2 * r / sigma ** 2
    delta_x = (b - a) / m
    delta_t = sigma ** 2 * T / (2 * nu_max)
    fidi_lambda = delta_t / delta_x ** 2
    lambda_tilde = (1 - 2 * fidi_lambda)

    ### range of underlying transformed stock prices
    x = np.arange(a, b + delta_x, delta_x)

    ### allocating memory
    w = np.zeros((m + 1, nu_max + 1))

    ### initial values equivalent to transformed payoff at maturity
    w[:, 0] = np.maximum(0, np.exp(x / 2 * (q + 1)) - np.exp(x / 2 * (q - 1)))

    ### loop over columns of matrix/time
    for nu in range(1, nu_max + 1):

        ### loop over rows/underlying (transformed) stock price
        ### note that we do not change the top and bottom row which are equal to zero all the time (simplified boundary conditions)
        for j in range(1, m):
            w[j, nu] = fidi_lambda * w[j - 1, nu - 1] + lambda_tilde * w[j, nu - 1] + fidi_lambda * w[j + 1, nu - 1]

        ### boundary condition for the next iteration
        #w[m, nu] = np.exp((q + 1) / 2 * b + (q + 1) ** 2 / 4 * nu * delta_t) - np.exp((q - 1) / 2 * b + (q - 1) ** 2 / 4 * nu * delta_t)

    ### retransfoming underlying stock prices
    S = K * np.exp(x[1:-1])

    ### transforming the solution of (5.1) into option prices
    V = K * w[1:-1, nu_max] * np.exp(-x[1:-1] / 2 * (q - 1) - sigma ** 2 / 2 * T * ((q - 1) ** 2 / 4 + q))
    return [S, V]


[S,V] = BS_EuCall_FiDi_Explicit(r, sigma, a, b, m, nu_max, T, K)


### BS-Formula
def EuCall_BlackScholes(t, S_t, r, sigma, T, K):
    d_1 = (math.log(S_t / K) + (r + 1 / 2 * math.pow(sigma, 2)) * (T - t)) / (sigma * math.sqrt(T - t))
    d_2 = d_1 - sigma * math.sqrt(T - t)
    Call = S_t * scipy.stats.norm.cdf(d_1) - K * math.exp(-r * (T - t)) * scipy.stats.norm.cdf(d_2)
    return Call

V_BS = np.zeros(len(S))

for i in range(0, len(S)):
    V_BS[i] = EuCall_BlackScholes(0, S[i], r, sigma, T, K)

plt.plot(S, V, label='Price with finite difference scheme')
plt.plot(S, V_BS, label='Price with BS-Formula')
plt.legend()
plt.show()

###### American Put Option


r = 0.05; sigma = 0.2; a = -0.7; b = 0.4; m = 100; nu_max = 2000; T = 1; K = 100

def BS_AmPut_FiDi_Explicit(r, sigma, a, b, m, nu_max, T, K):
    ### setting the parameters needed for the recursion
    q = 2 * r / sigma ** 2
    delta_x = (b - a) / m
    delta_t = sigma ** 2 * T / (2 * nu_max)
    fidi_lambda = delta_t / delta_x ** 2
    lambda_tilde = (1 - 2 * fidi_lambda)
    t = delta_t * np.arange(0, nu_max + 1)
    ### range of underlying transformed stock prices
    x = np.arange(a, b + delta_x, delta_x)

    ### allocating memory
    w = np.zeros((m + 1, nu_max + 1))
    #g = np.zeros((m + 1))
    
    ### initial values equivalent to transformed payoff at maturity CHANGE THIS FOR PUT 
    w[:, 0] = np.maximum(0, np.exp(x / 2 * (q - 1)) - np.exp(x / 2 * (q + 1)))

    ### loop over columns of matrix/time
    for nu in range(0, nu_max):
        for j in range(0, m):#CHANGE THIS FOR AMERICAN OPTIONS
            g = math.exp((q + 1) * (q + 1) * t[nu + 1] / 4) * np.maximum(np.exp(x * 0.5 * (q - 1))- np.exp(x * 0.5 * (q + 1)),np.zeros(m + 1))
            w[j, nu] = np.maximum(fidi_lambda * w[j - 1, nu - 1] + lambda_tilde * w[j, nu - 1] + fidi_lambda * w[j + 1, nu - 1], g[j])

        ### boundary condition for the next iteration
       # w[m, nu] = np.exp((q - 1) / 2 * b + (q - 1) ** 2 / 4 * t[nu] * delta_t) - np.exp((q + 1) / 2 * b + (q + 1) ** 2 / 4 * t[nu] * delta_t)

    ### retransfoming underlying stock prices
    S = K * np.exp(x[1:-1])

    ### transforming the solution of (5.1) into option prices
    #V = #K * w[1:-1, nu_max] * 
    V = np.exp(-x[1:-1] / 2 * (q - 1) - sigma ** 2 / 2 * T * ((q - 1) ** 2 / 4 + q))
    return S, V


[S, V3] = BS_AmPut_FiDi_Explicit(r, sigma, a, b, m, nu_max, T, K)

