% main.m

% Parameters
Pi0 = 1;            % Nominal value (100%)
T = 5;             % Maturity in years
delta_t = 1;        % Observation frequency in years
sigma = 0.3;        % Volatility
r = 0.02;           % Risk-free rate
Nmc = 100000;         % Number of simulations
timesteps = T / delta_t;  % Number of timesteps
S0_values = 80:1:120; % S0 varying from 80 to 120


% ---------- Fixed barriers
% Barriers
B_Ph_1 = 120;    % Fixed autocall barrier
B_Y_1 = 80;      % Fixed coupon barrier
B_Put_1 = 70;    % Fixed put barrier
B_Ph_2 = 100;    % Fixed autocall barrier
B_Y_2 = 70;      % Fixed coupon barrier
B_Put_2 = 60;    % Fixed put barrier
coeff_ph_1 = 0.1;
coeff_y_1 = 0.05;
coeff_ph_2 = 0.08;
coeff_y_2 = 0.05;
not_proportionnal = 0;
% Pricing cases
V_Case1_Put = Pricing_Phoenix(S0_values, Nmc, T, r, sigma, delta_t, Pi0, B_Ph_1, B_Y_1, B_Put_1, coeff_ph_1, coeff_y_1, 1, not_proportionnal);
V_Case1_Perf = Pricing_Phoenix(S0_values, Nmc, T, r, sigma, delta_t, Pi0, B_Ph_1, B_Y_1, B_Put_1, coeff_ph_1, coeff_y_1, 0, not_proportionnal);
V_Case2_Put = Pricing_Phoenix(S0_values, Nmc, T, r, sigma, delta_t, Pi0, B_Ph_2, B_Y_2, B_Put_2, coeff_ph_2, coeff_y_2, 1, not_proportionnal);
V_Case2_Perf = Pricing_Phoenix(S0_values, Nmc, T, r, sigma, delta_t, Pi0, B_Ph_2, B_Y_2, B_Put_2, coeff_ph_2, coeff_y_2, 0, not_proportionnal);

Graph_Pricing(S0_values, V_Case1_Put, V_Case1_Perf, V_Case2_Put, V_Case2_Perf);

% ---------- End fixed barriers

% ---------- Proportionnal barriers
% Barrier
B_Ph_coeff_1 = 1.2;   
B_Y_coeff_1 = 0.8;     
B_Put_coeff_1 = 0.7;  
B_Ph_coeff_2 = 1.0;
B_Y_coeff_2 = 0.7;
B_Put_coeff_2 = 0.6;

coeff_ph_1 = 0.1;
coeff_y_1 = 0.05;
coeff_ph_2 = 0.08;
coeff_y_2 = 0.05;
proportionnal = 1;
% Pricing cases
V_Case1_Put_Prop = Pricing_Phoenix(S0_values, Nmc, T, r, sigma, delta_t, Pi0, B_Ph_coeff_1, B_Y_coeff_1, B_Put_coeff_1, coeff_ph_1, coeff_y_1, 1, proportionnal);
V_Case1_Perf_Prop = Pricing_Phoenix(S0_values, Nmc, T, r, sigma, delta_t, Pi0, B_Ph_coeff_1, B_Y_coeff_1, B_Put_coeff_1, coeff_ph_1, coeff_y_1, 0, proportionnal);
V_Case2_Put_Prop = Pricing_Phoenix(S0_values, Nmc, T, r, sigma, delta_t, Pi0, B_Ph_coeff_2, B_Y_coeff_2, B_Put_coeff_2, coeff_ph_2, coeff_y_2, 1, proportionnal);
V_Case2_Perf_Prop = Pricing_Phoenix(S0_values, Nmc, T, r, sigma, delta_t, Pi0, B_Ph_coeff_2, B_Y_coeff_2, B_Put_coeff_2, coeff_ph_2, coeff_y_2, 0, proportionnal);

Graph_Pricing(S0_values, V_Case1_Put_Prop, V_Case1_Perf_Prop, V_Case2_Put_Prop, V_Case2_Perf_Prop);
% ---------- End proportional barriers

% Plot results


% Print specific prices for S0 = 100
S0_index = find(S0_values == 100);
fprintf('Prices for S0 = 100 with Fixed Barriers:\n');
fprintf('Case 1 - Put Payoff: %.6f\n', V_Case1_Put(S0_index));
fprintf('Case 1 - Performance Payoff: %.6f\n', V_Case1_Perf(S0_index));
fprintf('Case 2 - Put Payoff: %.6f\n', V_Case2_Put(S0_index));
fprintf('Case 2 - Performance Payoff: %.6f\n', V_Case2_Perf(S0_index));
