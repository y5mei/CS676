randn('state',100);
T = 0.50;
K = 102;
B = 100;
x = 15;
sigma = 0.20;
r = 0.03;
mu = 0.10;
S_init = 105;
N_sim = 100000; %number of simulations
delt = 0.5;     %temestep
N = T/delt;        %number of steps



drift = r-0.5*sigma*sigma;
sigma_sqrt_delt = sigma*sqrt(delt);

% initial values

S_old = zeros(N_sim, 1);
S_new = zeros(N_sim, 1);
payoff = zeros(N_sim, 1);

S_old(1:N_sim,1) = S_init;
payoff(1:N_sim,1) = x; 


for i = 1:N %  timestep loop
    
    % for each timestep, generate info for all simulations
    S_new(:,1) = S_old(:,1).*exp(drift*delt+sigma_sqrt_delt*randn(N_sim,1));
   
    z1 = log(S_new/K)/(sigma_sqrt_delt)+(sigma_sqrt_delt)/2;
    z2 = log(S_new/B)/(sigma_sqrt_delt)+(sigma_sqrt_delt)/2;  
    y1 = log((B*B/K)./S_new)/(sigma_sqrt_delt)+(sigma_sqrt_delt)/2; 
    y2 = log(B./S_new)/(sigma_sqrt_delt)+(sigma_sqrt_delt)/2;
    
    n1 = normcdf(-z1+sigma_sqrt_delt);
    n2 = normcdf(-z2+sigma_sqrt_delt);
    n3 = (S_new./B).*normcdf(y1-sigma_sqrt_delt);
    n4 = (S_new./B).*normcdf(y2-sigma_sqrt_delt);
    
    initial_price = x*exp(-r*T)*(n1-n2+n3-n4);
end

% replace any negative initial price with zero

greater_than_zero = (initial_price>0);
initial_price = initial_price.*greater_than_zero;
disp(mean(initial_price));

% The analytic solution when s0=105 is 0.0323



