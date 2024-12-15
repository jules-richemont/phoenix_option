% Pricing of a Down-and-In Barrier Put Option
Nmc = 1000;
S0 = 100;
K = 100;
r = 0.04;
sigma = 0.3;
T = 5;
B = 80;
N = 100;

price = Pricing_Barriere_Put_Down_In(Nmc, S0, K, r, sigma, B, T, N);
disp(['Price of a Down-and-In Barrier Put Option = ', num2str(price)]);

% Graph : 

S_min = 50;
S_max = 250;

Graph_Pricing_barriere_Put_Down_In(Nmc, K, r, sigma, B, T, S_min, S_max, N);


% Delta 

% Parameters
Nmc = 10000;
S0 = 100;
K = 100;
r = 0.05;
sigma = 0.2;
T = 1;
B = 90;
h = 0.1;

% Compute Delta
Delta = Compute_Delta(Nmc, S0, K, r, sigma, B, T, h);
% display
disp(['Delta = ', num2str(Delta)]);

% Graph Delta 
S_min = 50; 
S_max = 200; 

Graph_Delta(Nmc, K, r, sigma, B, T, S_min, S_max, h);
