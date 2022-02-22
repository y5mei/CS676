# This is the 3nd question sub question a for Assignment-2

# Compute Monte Carlo Cash-or-Nothing put barrier option
# European case
import math
import numpy
from scipy.interpolate import interp1d
import opstrat as op
import pandas as pd
from numpy.random import normal as randn
from numpy import mean
from scipy.stats import norm
import matplotlib.pyplot as plt
##########################################################
# Stock and Option Parameters
##########################################################


# Stock and Option Parameters
S0_list = list(range(100,122,2))
# S0 current stock price
T = 0.5   # T - expiry time
K = 102   # K - strike
B = 100   # B - down barrier
x = 15    # x - cash payout
M = 5000 # M - number of simulations
# N = 10   # N - number of timesteps
r = 0.03  # r - interest rate
sigma = 0.2  # sigma - volatility
# delt = T/N  # timestep


def Sn_plus_one(Sn: float, random_num: float, delt: float = T) -> float:
    # Using BS, calculate Sn+1 from Sn
    # Sn is the current Stock Price
    # random_num is a random number draw from a standard normal distribution
    # delt is the timestep, default equals to T
    # Sn_plus_one(S0) returns ST
    # Sn_plus_one(Sn) returns Sn+1
    drift = (r - sigma * sigma / 2) * delt
    sigma_sqrt_delt = sigma * math.sqrt(delt)
    return Sn*math.exp(drift+sigma_sqrt_delt*random_num)


def V0(St: float) -> float:
    # calculate V(S,0)
    z1 = math.log(St/K)/(sigma*math.sqrt(T))+0.5*sigma*math.sqrt(T)
    z2 = math.log(St/B)/(sigma*math.sqrt(T))+0.5*sigma*math.sqrt(T)
    y1 = math.log((B*B)/(St*K))/(sigma*math.sqrt(T))+0.5*sigma*math.sqrt(T)
    y2 = math.log(B/St)/(sigma * math.sqrt(T)) + 0.5 * sigma * math.sqrt(T)

    x1 = -z1+sigma*math.sqrt(T)
    x2 = -z2+sigma*math.sqrt(T)
    x3 = y1-sigma*math.sqrt(T)
    x4 = y2-sigma*math.sqrt(T)

    return x*math.exp(-r*T)*(norm.cdf(x1)-norm.cdf(x2)+(St/B)*norm.cdf(x3)-(St/B)*norm.cdf(x4))



def Exacl_Solution(S0_list):
    initial_option_value = []
    final_stock_price = []
    for s0 in S0_list:
        s0_list = [s0]*M
        rand_list = randn(0,1,M)
        st_list = list(map(Sn_plus_one, s0_list, rand_list)) # list of simulated final stock price
        final_stock_price.append(mean(st_list))                             # average of simulated final stock price
        v0_list = list(map(V0, st_list))                     # calculated initial option price
        v0 = mean(v0_list)
        initial_option_value.append(v0)

    return initial_option_value

initial_option_value = Exacl_Solution(S0_list)

## make the plot of option fair value against inital stock price
fig,ax = plt.subplots(1)
ax.plot(S0_list,initial_option_value, color='green', marker='o', linestyle='dashed',
     linewidth=2, markersize=10)

font = {'family': 'serif',
        'color':  'darkred',
        'weight': 'normal',
        'size': 16,
        }

ax.set_xlabel("Initial Stock Price",fontdict=font)
ax.set_ylabel("Initial Option Fair Value",fontdict = font)
plt.show()


FIXED_S0 = 105
exact_solution = Exacl_Solution([FIXED_S0])


print("The Exact Solution for S0 = 105 is {0:.5f}".format(exact_solution[-1]))



