import numpy as np 


def Delta(beta, alpha, d):
    return_value = d*(np.sin(beta)+np.sin(alpha))
    #if return_value < 0.00000001 and return_value > -0.00000001:
    #    return_value = 0
    #print('Delta, ',return_value)
    return return_value

def delta(alpha, sigma, beta):
    return_value = (np.sin(alpha-sigma)+np.sin(beta-sigma))/np.cos(sigma)
    #if return_value < 0.00000001 and return_value > -0.00000001:
    #    return_value = 0
    #print('delta, ', return_value)
    return return_value

def T_fct(lambda_,  N, d, alpha, beta, sigma):
    return  (np.sinc((np.pi*d*delta(alpha, sigma, beta))/lambda_ ) *((np.sin(np.pi*N*Delta(beta, alpha, d)/lambda_)) /(N*np.sin(np.pi*Delta(beta, alpha, d)/lambda_))))**2

    
def func_beta(m, lambda_, d, alpha):
    x = m*lambda_/d -np.sin(alpha)
    for i in range(len(x)):
        if x[i] <= 1 and x[i] >= -1:
            x[i]=np.arcsin(x[i]) 
        else:
            x[i]=np.NaN
    return x

def slitsize(d, lambda_, alpha, f):
    return f*0.1*10**(-9)/(d*np.cos(func_beta(1,lambda_, d, alpha)))