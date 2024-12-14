% main.m

% Parameters
Pi0 = 1;            % Nominal value (100%)
T = 10;             % Maturity in years
delta_t = 1;        % Observation frequency in years
sigma = 0.3;        % Volatility
r = 0.02;           % Risk-free rate
N = 100000;         % Number of simulations
timesteps = T / delta_t;
S0_values = 80:1:120; % S0 varying from 80 to 120

% Barriers
B_Ph = 120;    % Fixed autocall barrier
B_Y = 80;      % Fixed coupon barrier
B_Put = 70;    % Fixed put barrier
coeff_ph_1 = 0.1;
coeff_y_1 = 0.05;
coeff_ph_2 = 0.08;
coeff_y_2 = 0.05;

% Pricing cases
V_Case1_Put = Pricing_Phoenix(S0_values, N, T, r, sigma, delta_t, Pi0, B_Ph, B_Y, B_Put, coeff_ph_1, coeff_y_1, 1);
V_Case1_Perf = Pricing_Phoenix(S0_values, N, T, r, sigma, delta_t, Pi0, B_Ph, B_Y, B_Put, coeff_ph_1, coeff_y_1, 0);
V_Case2_Put = Pricing_Phoenix(S0_values, N, T, r, sigma, delta_t, Pi0, B_Ph, B_Y, B_Put, coeff_ph_2, coeff_y_2, 1);
V_Case2_Perf = Pricing_Phoenix(S0_values, N, T, r, sigma, delta_t, Pi0, B_Ph, B_Y, B_Put, coeff_ph_2, coeff_y_2, 0);


% Plot results
Graph_Pricing(S0_values, V_Case1_Put, V_Case1_Perf, V_Case2_Put, V_Case2_Perf);


% Print specific prices for S0 = 100
S0_index = find(S0_values == 100);
fprintf('Prices for S0 = 100 with Fixed Barriers:\n');
fprintf('Case 1 - Put Payoff: %.6f\n', V_Case1_Put(S0_index));
fprintf('Case 1 - Performance Payoff: %.6f\n', V_Case1_Perf(S0_index));
fprintf('Case 2 - Put Payoff: %.6f\n', V_Case2_Put(S0_index));
fprintf('Case 2 - Performance Payoff: %.6f\n', V_Case2_Perf(S0_index));
