"""draft"""
#%%
import math
import numpy as np
from scipy.stats import norm

def BSM_analytical(type, S0, X, T, r, sigma, div):
    # zero strike call = price of underlying without dividends
    if X==0:
        bsm = S0*np.exp(-div*T)
    else:
        d_1 = (np.log(S0/X)+(r+sigma**2 /2)*T)/(sigma*np.sqrt(T))
        d_2 = (np.log(S0/X)+(r-sigma**2 /2)*T)/(sigma*np.sqrt(T))
    
    if type=='Call':
        bsm = S0*norm.cdf(d_1)-X*np.exp(-r*T)*norm.cdf(d_2)
    elif type=='Put':
        bsm = X*np.exp(-r*T)*norm.cdf(-d_2) - S0*norm.cdf(-d_1)
    
    else:
        print('Type not valid')
    
    return bsm
 
 
def guaranteeCertificate(S, X, T, r, sigma, nominal, participation, ratio):
    
    S0 = S
    # zero coupon bond, nominal corresponds to guaranteed payoff level
    zcb = nominal * np.exp(-r * T)
    
    if T==0:
        long_call = max(S - X, 0) * participation
    else:
        long_call = BSM_analytical('Call', S0, X, T, r, sigma, 0) * participation
        
    value = (zcb + long_call) * ratio
    
    return value

#test
test_price = guaranteeCertificate(110, 90, 1, 0.045, 0.4, 110, 0.8, 1)
test_price
# %%
