## Hedging Finite Difference Explicit 


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
        #w[m, nu] = np.exp((q + 1) / 2 * b + (q + 1) ** 2 / 4 * nu * delta_t) - np.exp(
         #   (q - 1) / 2 * b + (q - 1) ** 2 / 4 * nu * delta_t)

    ### retransfoming underlying stock prices
    S = K * np.exp(x[1:-1])

    ### transforming the solution of (5.1) into option prices
    V = K * w[1:-1, nu_max] * np.exp(-x[1:-1] / 2 * (q - 1) - sigma ** 2 / 2 * T * ((q - 1) ** 2 / 4 + q))

    ### using finite difference quotient to caclcute the phi_1(0), i.e. the first order derivative of the option price V(0) w.r.t. to S(0)
    phi = np.diff(V)/ np.diff(S)
    return [phi, S, V]


[phi,S, V] = BS_EuCall_FiDi_Explicit(r, sigma, a, b, m, nu_max, T, K)


### Exact hedge in the BS-model
def EuCall_BlackScholes_hedge(t, S_t, r, sigma, T, K):
    d_1 = (math.log(S_t / K) + (r + 1 / 2 * math.pow(sigma, 2)) * (T - t)) / (sigma * math.sqrt(T - t))
    return scipy.stats.norm.cdf(d_1)

V_BS_hedge = np.zeros(len(S)-1)
### applying the BS-Formula to each underlying stock price
### note that we do set the stock prices only indirectly through the parameters a,b and m
for i in range(0, (len(S)-1)):
    V_BS_hedge[i] = EuCall_BlackScholes_hedge(0, S[i], r, sigma, T, K)

plt.plot(S[0:-1], phi, label='Price with finite difference scheme')
plt.plot(S[0:-1], V_BS_hedge, label='Price with BS-Formula')
plt.legend()
plt.show()
