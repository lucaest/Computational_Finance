
""" 
implementation as in og paper
"""
#%%

import math
import numpy as np
np.random.seed(93)

S0 = 1.  
K = 1.1  
T = 3.0  
r = 0.06 
sigma = 0.2 

I = 8
M = 4
dt = T / (M-1)
df = math.exp(-r * dt)

def payoff(type, S, K):
    if type=='Call':
        value = np.maximum(S - K, 0)
    elif type=='Put':
        value = np.maximum(K - S, 0)
    else:
        print('Type not valid')

    return value

# for general values 
S = np.zeros((I, M))
S[:, 0] = S0
X = np.random.normal(0, 1, (I, M))

for i in range(0, I):
    for j in range(1, M):
        S[i, j] = S0 * math.exp((r-sigma**2/2)*dt + sigma*math.sqrt(dt)*X[i, j])

# for paper values, overwrite
S = np.array([
           [1.00, 1.09, 1.08, 1.34],
           [1.00, 1.16, 1.26, 1.54],
           [1.00, 1.22, 1.07, 1.03],
           [1.00, 0.93, 0.97, 0.92],
           [1.00, 1.11, 1.56, 1.52],
           [1.00, 0.76, 0.77, 0.90],
           [1.00, 0.92, 0.84, 1.01],
           [1.00, 0.88, 1.22, 1.34]])
        
H = payoff('Put', S, K)      
V = np.zeros((I, M))
V[:,-1] = H[:,-1]

for t in range(M-2, 0, -1):
    rg = np.polyfit( S[:, t], V[:, t+1] * df, 2)  
    C = np.polyval( rg, S[:,t] )        
    V[:, t] = np.where(H[:, t] > C, H[:, t], V[:, t+1] * df) 
    
V0 = np.mean(V[:, 1]) * df 
print('option value: ' + str(V0))
# check: 0.1144
# %%
