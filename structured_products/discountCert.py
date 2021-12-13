"""draft"""
#%%
import math
import numpy as np
from scipy.stats import norm

def BSM_analytical(type, S0, X, T, r, sigma):
    # FIX division by zero 
    if X==0:
        d_1 = (np.log(S0)+(r+sigma**2 /2)*T)/(sigma*np.sqrt(T))
        d_2 = (np.log(S0)+(r-sigma**2 /2)*T)/(sigma*np.sqrt(T))
    
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

def discountCertificate(S, X, T, r, sigma):
    
    # zero strike call position
    zero_strike_call = BSM_analytical(type='Call', S0=S, X=0, T=T, r=r, sigma=sigma)
    
    # short call position
    short_call = BSM_analytical(type='Call', S0=S, X=X, T=T, r=r, sigma=sigma)
    
    value = zero_strike_call - short_call
    
    return value



# test run:
# test_price = discountCertificate(100, 80, 9/12, 0.06, 0.35)
# test_price

# %%
