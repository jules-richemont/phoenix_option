% Pricing of a Down-and-In Barrier Put Option
Nmc = 1000;
S0 = 100;
K = 100;
r = 0.05;
sigma = 0.2;
T = 1;
B = 90;

price = Pricing_Barriere_Put_Down_In(Nmc, S0, K, r, sigma, B, T);


% Graph : 

Nmc = 10000;
K = 100;
r = 0.05;
sigma = 0.2;
B = 90;
T = 1;
S_min = 50;
S_max = 200;

Graph_Pricing_barriere_Put_Down_In(Nmc, K, r, sigma, B, T, S_min, S_max);


% Delta 

% Parameters
Nmc = 10000;
S0 = 100; % Current asset price
K = 100; % Strike price
r = 0.05; % Risk-free rate
sigma = 0.2; % Volatility
T = 1; % Maturity
B = 90; % Barrier level
h = 0.1; % Step size for Delta computation

% Compute Delta
Delta = Compute_Delta(Nmc, S0, K, r, sigma, B, T, h);
% display
disp(['Delta = ', num2str(Delta)]);

% Graph Delta 
S_min = 50; % Minimum S0 for the plot
S_max = 200; % Maximum S0 for the plot

Graph_Delta(Nmc, K, r, sigma, B, T, S_min, S_max, h);
