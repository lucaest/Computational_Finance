"""draft"""
#%%
import math
import numpy as np
from scipy.stats import norm

def BSM_analytical(type, S0, X, T, r, sigma, div):
    # zero strike call = price of underlying without dividends
    if type=='ZeroStrike':
        bsm = S0*np.exp(-r*T)
    elif type=='Call':
        d_1 = (np.log(S0/X)+(r+sigma**2 /2)*T)/(sigma*np.sqrt(T))
        d_2 = (np.log(S0/X)+(r-sigma**2 /2)*T)/(sigma*np.sqrt(T))
        bsm = S0*norm.cdf(d_1)-X*np.exp(-r*T)*norm.cdf(d_2)
    elif type=='Put':
        d_1 = (np.log(S0/X)+(r+sigma**2 /2)*T)/(sigma*np.sqrt(T))
        d_2 = (np.log(S0/X)+(r-sigma**2 /2)*T)/(sigma*np.sqrt(T))
        bsm = X*np.exp(-r*T)*norm.cdf(-d_2) - S0*norm.cdf(-d_1)
    else:
        print('Type not valid')
    
    return bsm


class barrierOptions():
    """ 
    X: strike
    H: barrier
    K: cash rebate
    b: annual cost of carry
    """
    def standardBarrierOption(type, barrierType, S, X, H, K, T, r, b, sigma, div):
    
        if type == 'Call':
            if barrierType == 'downAndIn' or barrierType == 'downAndOut':
                phi = 1
                eta = 1
            elif barrierType == 'UpAndIn' or barrierType == 'UpAndOut':
                phi = 1
                eta = -1
            else:
                print('invalid barrier type')
        elif type == 'Put':
            if barrierType == 'downAndIn' or barrierType == 'downAndOut':
                phi = -1
                eta = 1
            elif barrierType == 'UpAndIn' or barrierType == 'UpAndOut':
                phi = -1
                eta = -1
            else:
                print('invalid barrier type')        
        else:
            print('invalid option type')
        
        mu = (b - sigma**2 / 2) / sigma**2
        lbd = np.sqrt(mu**2 + 2 * r / sigma**2)
        
        x1 = math.log(S / X) / (sigma * np.sqrt(T)) + (1 + mu) * sigma * np.sqrt(T)
        x2 = math.log(S / H) / (sigma * np.sqrt(T)) + (1 + mu) * sigma * np.sqrt(T)
        
        y1 = math.log(H**2 / (S * X)) / (sigma * np.sqrt(T)) + (1 + mu) * sigma * np.sqrt(T)
        y2 = math.log(H / S) / (sigma * np.sqrt(T)) + (1 + mu) * sigma * np.sqrt(T)
        
        z  = math.log(H / S) / (sigma * np.sqrt(T)) + lbd * sigma * np.sqrt(T)

        A = (phi * S * math.exp((b - r) * T) * norm.cdf(phi * x1) - phi * X * math.exp(-r * T) * norm.cdf(phi * x1 - phi * sigma * np.sqrt(T)))
        B = (phi * S * math.exp((b - r) * T) * norm.cdf(phi * x2) - phi * X * math.exp(-r * T) * norm.cdf(phi * x2 - phi * sigma * np.sqrt(T)))
        
        C = (phi * S * math.exp((b - r) * T) * (H / S)**(2 * (mu + 1)) * norm.cdf(eta * y1) - phi * X * math.exp(-r * T) * (H / S)**(2 * mu) * norm.cdf(eta * y1 - eta * sigma * np.sqrt(T)))
        D = (phi * S * math.exp((b - r) * T) * (H / S)**(2 * (mu + 1)) * norm.cdf(eta * y2) - phi * X * math.exp(-r * T) * (H / S)**(2 * mu) * norm.cdf(eta * y2 - eta * sigma * np.sqrt(T)))
        
        E = K * math.exp(-r * T) * (norm.cdf(eta * x2 - eta * sigma * np.sqrt(T)) - (H / S)**(2 * mu) * norm.cdf(eta * y2 - eta * sigma * np.sqrt(T)))
        F = K * (H / S)**(mu + lbd) * norm.cdf(eta * z) + (H / S)**(mu - lbd) * norm.cdf(eta * z - 2 * eta * lbd * sigma * np.sqrt(T))       
        
        if type == 'Call':
            if barrierType == 'downAndIn':
                if X > H:
                    standardBarrier = C + E
                elif X < H:
                    standardBarrier = A - B + D + E
            elif barrierType == 'UpAndIn':
                if X > H:
                    standardBarrier = A + E
                elif X < H:
                    standardBarrier = B - C + D + E
            elif barrierType == 'downAndOut':
                if X > H:
                    standardBarrier = A - C + F
                elif X < H:
                    standardBarrier = B - D + F
            elif barrierType == 'UpAndOut':
                if X > H:
                    standardBarrier = F
                elif X < H:
                    standardBarrier = A - B + C - D + F
            else:
                print('invalid barrier type')
        elif type == 'Put':
            if barrierType == 'downAndIn':
                if X > H:
                    standardBarrier = B - C + D + E
                elif X < H:
                    standardBarrier = A + E
            elif barrierType == 'UpAndIn':
                if X > H:
                    standardBarrier = A - B + D + E
                elif X < H:
                    standardBarrier = C + E
            if barrierType == 'downAndOut' :
                if X > H:
                    standardBarrier = A - B + C - D + F
                elif X < H:
                    standardBarrier = F
            elif barrierType == 'UpAndOut':
                if X > H:
                    standardBarrier = B - D + F
                elif X < H:
                    standardBarrier = A - C + F
            else:
                print('invalid barrier type')        
        else:
            print('invalid option type')
          
        return standardBarrier
    

def bonusCertificate(S, X, H, T, r, sigma, ratio):
    
    # zero strike call position
    zero_strike_call = BSM_analytical('ZeroStrike', S, 0, T, r, sigma, 0)
    
    # down and out put option
    barrierPut = barrierOptions.standardBarrierOption('Put', 'downAndOut', S, X, H, 0, T, r, 0, sigma, 0)
    
    value = (zero_strike_call - barrierPut) * ratio
    
    return value

#test
test_price = bonusCertificate(50, 60, 35, 2, 0.02, 0.14, 1)
test_price
