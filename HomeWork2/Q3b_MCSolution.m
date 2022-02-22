function result = fast()
delt_list = [200, 400, 800, 1600, 3200, 6400];
M_list = 1000:3000:100000;

len = size(M_list,2);

result = [];
% mylengend = cell(1, len);
% i = 1;
% line_color = ['b' 'g' 'y' 'c' 'm' 'r'];

for d = delt_list
    d_list = zeros(len,1);
    d_list(1:len,1) = 1/d;
    d_list = d_list';
    
    disp(d_list);
    disp(M_list);
    
%     mylengend{i} = sprintf("1/%d", d);
%     i = i+1;
    
    r = arrayfun(@(x,y) MC_Solution(x,y), d_list, M_list);
%     disp(r);
    result  = [result; r];
end

% Exact Solution
exact_solution = zeros(len,1);
exact_solution(1:len,1) = 0.0322644;
exact_solution = exact_solution';
% result  = [result; exact_solution];

h = plot(M_list,result,'--o',M_list,exact_solution,'LineWidth',3);
xtickangle(45);
legend('1/200','1/400','1/800','1/1600','1/3200','1/6400','Exact');
end


function initial_value = MC_Solution(delt,N_sim)

delt = delt;     %temestep
N_sim = N_sim; %number of simulations
randn('state',100);
T = 0.50;
K = 102;
B = 100;
x = 15;
sigma = 0.20;
r = 0.03;
mu = 0.10;
S_init = 105;
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
    
    % make an arry, not_hits_barrier,
    %     1) if any new_value hits B, payoff will be zero, we mark the
    %     elements in arry as 0;
    %     2) if any new_value not hits B, payoff will not be zero, we mark the
    %     elements in array as 1; 
    not_hits_barrier = (S_new > B);
    
    % multiply payoff with the boolean array, every elements was zero
    % before, will always be zero;
    payoff = payoff.* not_hits_barrier;
%     S_new = S_new.*not_hist_barrier;
    
    S_old(:,1) = S_new(:,1); % end of generating all data for all simulations for this timestep
end


% now, check the stock value at T
sk_smaller_than_k = (S_old < K);
payoff = payoff.* sk_smaller_than_k;

initial_value = exp(-r*T)*mean(payoff);

% disp(initial_value)

end

