# This is the 1st quesiton for Assignment-2

# Compute Black-Scholes option value suing a binomial tree
# European case
import math
import numpy
from scipy.interpolate import interp1d
import opstrat as op
import pandas as pd

##########################################################
# Stock and Option Parameters
##########################################################

# Stock and Option Parameters
S0 = 10  # S0 current stock price
K = 10  # K - strike
T = 1.0  # T - expiry time
r = 0.03  # r - interest rate
sigma = 0.2  # sigma - volatility
opttype = 0  # opttype - 0 for a call, otherwise a put
delta = 0.00625  # time_step of the lattice


# Dividend parameters
D0 = 0.5  # D0 the floor of dividend
rho = 0  # constant dividend rate
dividendtype = 0  # dividendtype - 0 for no dividend, 1 for regular dividend at T/3
td = T / 3  # dividend payment time


##########################################################
# helper functions to
# 1) Calculate dividend payment,
# 2) Calculate the ex_stock price (S_d_t^+)
##########################################################

# Calculate the amount of dividend payment D
def D(stock_price_before_dividend: float, rho: float) -> float:
    return max(rho * stock_price_before_dividend, D0)


# Calculate the stock value ex dividend (after pay the dividend)
def ex_dividend_stock_value(stock_price_before_dividend: float, rho: float) -> float:
    return stock_price_before_dividend - D(stock_price_before_dividend, rho)

##########################################################
# The main method
#
##########################################################

def binomial_lattice_pricing(opttype:int = 0, dividendtype:int = 0, rho:float = 0, delta:float = 0.05) -> float:

    # Calculate the size of lattice and the node of dividend
    Nsteps = int(T / delta)  # Nsteps - number of timesteps
    NDsteps = round(td / delta)  # NDsteps - the lattice node where the dividend get paid

    # Tree parameters
    u = math.exp(sigma * math.sqrt(delta))
    d = 1 / u
    a = math.exp(r * delta)
    p = (a - d) / (u - d)

    # Payoff at t=T
    W = [S0 * math.pow(u, i) * math.pow(d, Nsteps - i) for i in range(Nsteps + 1)]

    if opttype == 0:
        # if this is a call option
        # print("This is a call option")
        W = list(map(lambda s: max(s - K, 0), W))
    else:
        # if this is a put option
        # print("This is a put option")
        W = list(map(lambda s: max(K - s, 0), W))

    # Backward recursion

    for n in range(Nsteps - 1, -1, -1):
        # print("This is month ", n)
        for j in range(n + 1):
            W[j] = (1 / a) * (p * W[j + 1] + (1 - p) * W[j])
        if n == NDsteps and dividendtype:
            # print(
            #     "There are totally {:d} types in this lattice, and the dividend paid at step {:d}".format(Nsteps,
            #                                                                                               NDsteps))
            # calculate the stock price before dividend
            S = [S0 * math.pow(u, i) * math.pow(d, NDsteps - i) for i in range(NDsteps + 1)]
            # calculate the stock price ex paying the dividend
            list_rho = [rho]*len(S)
            S_ex = list(map(ex_dividend_stock_value, S, list_rho))
            # Build the interpolation function:
            interp_func = interp1d(S, W[:NDsteps + 1], fill_value="extrapolate")
            # calculate Option Value at t^-
            W = interp_func(S_ex)
    return W[0]





print("===============================================================================================")
print("==========================          Assignment Q1-(a)     =====================================")
print("===============================================================================================")

# BS Model
# Declare parameters
# K = 10  # spot price
# St = 10  # current stock price
# r = 3  # 3% risk free rate
# t = 365  # time to expiry, 365 days = 1 year
# v = 20  # volatility
# type = 'c'  # Option type call
print("========================== Standard Result  from blsprice ======================================")
bsm_call = op.black_scholes(K=10, St=10, r=3, t=365, v=20, type='c')
print("blsprice call option value is:", bsm_call['value']["option value"])  # The standard blsprice call value should be 0.9413403383853005
bsm_put = op.black_scholes(K=10, St=10, r=3, t=365, v=20, type='p')
print("blsprice put option value is:",bsm_put['value']["option value"])  # The standard blsprice put value should be 0.6457956738703823

print("========================== My Calculation for Call Option ======================================")

delta_t_list = []
value_list = []
change_list = ["NA"]
ratio_list = ["NA", "NA"]
nsteps_list = []
error_list = []
call_blsprice_value = 0.9413403383853005

delta_t = 0.05
for i in range(8):
    delta_t_list.append(delta_t)
    value_list.append(binomial_lattice_pricing(opttype=0, dividendtype=0, delta=delta_t)) # call
    nsteps_list.append(int(1/delta_t))
    error_list.append((abs(value_list[-1] - call_blsprice_value) / call_blsprice_value) * 100)
    if i > 0:
        change_list.append(value_list[i] - value_list[i - 1])
    if i > 1:
        ratio_list.append(change_list[i-1]/change_list[i])
    delta_t/=2

# dictionary of lists
call_option_dict = {'Delta_t': delta_t_list, 'Value': value_list, 'Change': change_list, "Ratio":ratio_list, "Nsteps":nsteps_list, "Error%": error_list}
df = pd.DataFrame(call_option_dict)
print(df)
df.to_csv("call_option_no_dividend.csv")
print("===============================================================================================")
print("========================== My Calculation for Put Option ======================================")

delta_t_list = []
value_list = []
change_list = ["NA"]
ratio_list = ["NA", "NA"]
nsteps_list = []
error_list = []
put_blsprice_value = 0.6457956738703823

delta_t = 0.05
for i in range(8):
    delta_t_list.append(delta_t)
    value_list.append(binomial_lattice_pricing(opttype=1, dividendtype=0, delta=delta_t)) # put
    nsteps_list.append(int(1 / delta_t))
    error_list.append ((abs(value_list[-1]-put_blsprice_value)/put_blsprice_value)*100)
    if i > 0:
        change_list.append(value_list[i] - value_list[i - 1])
    if i > 1:
        ratio_list.append(change_list[i-1]/change_list[i])
    delta_t/=2

# dictionary of lists
put_option_dict = {'Delta_t': delta_t_list, 'Value': value_list, 'Change': change_list, "Ratio":ratio_list,"Nsteps":nsteps_list, "Error%": error_list}
df = pd.DataFrame(put_option_dict)
print(df)
df.to_csv("put_option_no_dividend.csv")
print("===============================================================================================")

print("===============================================================================================")
print("==========================          Assignment Q1-(b)     =====================================")
print("===============================================================================================")
rho_list = [0, 0.02, 0.04, 0.08]
call_option_values_list = []
put_option_values_list = []
for rho in rho_list:
    call_option_values_list.append(binomial_lattice_pricing(opttype=0, dividendtype=1, delta=0.005, rho = rho))
    put_option_values_list.append(binomial_lattice_pricing(opttype=1, dividendtype=1, delta=0.005, rho=rho))

question2_result_dict = {'Rho': rho_list, "Call_Value": call_option_values_list, "Put_Value": put_option_values_list}
df = pd.DataFrame(question2_result_dict)
print(df)