"""draft"""
#%%
import math
import numpy as np
from scipy.stats import norm


def BSM_analytical(type, S0, X, T, r, sigma):
    d_1 = (np.log(S0/X)+(r+sigma**2 /2)*T)/(sigma*np.sqrt(T))
    d_2 = (np.log(S0/X)+(r-sigma**2 /2)*T)/(sigma*np.sqrt(T))
    
    if type=='Call':
        bsm = S0*norm.cdf(d_1)-X*np.exp(-r*T)*norm.cdf(d_2)
    elif type=='Put':
        bsm = X*np.exp(-r*T)*norm.cdf(-d_2) - S0*norm.cdf(-d_1)
    
    else:
        print('Type not valid')
    
    return bsm


def cashornothingOption(type, S, X, K, b, T, r, sigma):
    # K is the amount that is paid if option is ITM, else 0
    S0=S
    # get d_2 from BSM_analytical
    d_2 = (np.log(S0/X) + (b-sigma**2 /2)*T)/(sigma*np.sqrt(T))

    if type=='Call':
        cash_or_nothing = K * np.exp(-r*T) * norm.cdf(d_2)
    
    elif type=='Put':
        cash_or_nothing = K * np.exp(-r*T) * norm.cdf(-d_2)
        
    else:
        print('Type not valid')
    
    
    return cash_or_nothing

def expressCertificate(S, S0, X, T, r, sigma):
    
    # long zero bond position
    zero_bond = S0 * np.exp(-r*T)
    
    # short cash or nothing put position
    cash_or_nothing = cashornothingOption(type='Put', S=S, X=X, K=S0-X, b=0, T=T, r=r, sigma=sigma)
    
    # short put position
    short_put = BSM_analytical(type='Put', S0=S0, X=X, T=T, r=r, sigma=sigma)
    
    value = zero_bond - cash_or_nothing - short_put
    
    return value

# %%
# test run:
# test_price = expressCertificate(80, 100, 70, 1, 0.045, 0.4)
# test_price
# %%
