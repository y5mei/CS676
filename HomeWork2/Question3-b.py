# This is the 3nd question sub question b for Assignment-2

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
M = 30000 # M - number of simulations
r = 0.03  # r - interest rate
sigma = 0.2  # sigma - volatility
# delt = T/N  # timestep
# N = T/delt   # N - number of timesteps

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



def MC_Solution(delt, M, s0 = 105):
    # calculate the initial option price using a time step delt, and a limit of simulation num. M
    s_list = [s0] * M
    delt_list = [delt] * M
    payoff_list = [x]*M

    Nsteps = int(T/delt)
    for n in range(Nsteps+1):
        rand_list = randn(0, 1, M)
        s_list = list(map(Sn_plus_one, s_list, rand_list, delt_list)) # list of simulated final stock price
        for i in range(M):
            if s_list[i] <= B:
                payoff_list[i] = 0 # payoff = 0 if it hits B

    # check the final stock price with strike price
    for i in range(M):
        if s_list[i] >= K:
            payoff_list[i] = 0 # payoff = 0 if stock price higher than K


    # calculate the average payoff
    average_pay_off = mean(payoff_list)
    initial_option_value = math.exp(-r*T)*average_pay_off
    return initial_option_value

# initial_option_value = MC_Solution(105, 1 / 200, M)

M_LIST = list(range(1000, 103000, 6000))
DELTA_T_LIST = [1/200, 1/400, 1/800, 1/1600, 1/3200, 1/6400]
# DELTA_T_LIST = [1/200, 1/400]


from datetime import datetime

now = datetime.now()
current_time = now.strftime("%H:%M:%S")
print("Start Time =", current_time)

result = []
for delt in DELTA_T_LIST:
    delt_list = [delt]*len(M_LIST)
    mc_solution = list(map(MC_Solution,delt_list, M_LIST))
    result.append(mc_solution)



now = datetime.now()
current_time = now.strftime("%H:%M:%S")
print("End Time =", current_time)

## make the plot of option fair value against inital stock price
## and save result to a CSV
fig,ax = plt.subplots(1)
put_option_dict = {}
for i in range(len(result)):
    line = result[i]
    label = "delta_t = "+str(DELTA_T_LIST[i])
    put_option_dict[label] = line # add to dataframe
    ax.plot(M_LIST, line, linestyle='dashed', marker='o', markersize=5, label=label) # make the plot


df = pd.DataFrame(put_option_dict)
print(df)
df.to_csv("Q3_put_option_MC_Method.csv")

font = {'family': 'serif',
        'color':  'darkred',
        'weight': 'normal',
        'size': 16,
        }

ax.set_xlabel("Number of Simulation",fontdict=font)
ax.set_ylabel("Initial Option Fair Value",fontdict = font)
plt.legend()
plt.show()







