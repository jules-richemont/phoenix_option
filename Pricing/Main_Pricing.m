% main.m

% Parameters
Pi0 = 1;            
T = 10;             
delta_t = 1;        
sigma = 0.3;        
r = 0.02;           
Nmc = 10000;    
N = 100;            
S0_values = 0:0.5:150; 


% ---------- Fixed barriers
% Barriers
B_Ph_1 = 120; % Autocall barrier
B_Y_1 = 80;% Coupon barrier
B_Put_1 = 70;% Put barrier   
B_Ph_2 = 100;%     
B_Y_2 = 70;      
B_Put_2 = 60;    
coeff_ph_1 = 0.1;
coeff_y_1 = 0.05;
coeff_ph_2 = 0.08;
coeff_y_2 = 0.05;
not_proportionnal = 0;
% Pricing cases
V_Case1_Put = Pricing_Phoenix(S0_values, Nmc, T, r, sigma, delta_t, Pi0, B_Ph_1, B_Y_1, B_Put_1, coeff_ph_1, coeff_y_1, 1, not_proportionnal, N);
V_Case1_Perf = Pricing_Phoenix(S0_values, Nmc, T, r, sigma, delta_t, Pi0, B_Ph_1, B_Y_1, B_Put_1, coeff_ph_1, coeff_y_1, 0, not_proportionnal, N);
V_Case2_Put = Pricing_Phoenix(S0_values, Nmc, T, r, sigma, delta_t, Pi0, B_Ph_2, B_Y_2, B_Put_2, coeff_ph_2, coeff_y_2, 1, not_proportionnal, N);
V_Case2_Perf = Pricing_Phoenix(S0_values, Nmc, T, r, sigma, delta_t, Pi0, B_Ph_2, B_Y_2, B_Put_2, coeff_ph_2, coeff_y_2, 0, not_proportionnal, N);
% Print specific prices for S0 = 100
S0_index = find(S0_values == 100);
fprintf('Prices for S0 = 100 with Fixed Barriers:\n');
fprintf('Case 1 - Put Payoff: %.6f\n', V_Case1_Put(S0_index));
fprintf('Case 1 - Performance Payoff: %.6f\n', V_Case1_Perf(S0_index));
fprintf('Case 2 - Put Payoff: %.6f\n', V_Case2_Put(S0_index));
fprintf('Case 2 - Performance Payoff: %.6f\n', V_Case2_Perf(S0_index));
Graph_Pricing(S0_values, V_Case1_Put, V_Case1_Perf, V_Case2_Put, V_Case2_Perf, not_proportionnal);

% ---------- End fixed barriers

% ---------- Proportionnal barriers

% Parameters
Pi0 = 1;            
T = 10;             
delta_t = 1;       
sigma = 0.3;    
r = 0.02;        
Nmc = 10000;        
N = 100;            
S0_values = 80:1:120;
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
V_Case1_Put_Prop = Pricing_Phoenix(S0_values, Nmc, T, r, sigma, delta_t, Pi0, B_Ph_coeff_1, B_Y_coeff_1, B_Put_coeff_1, coeff_ph_1, coeff_y_1, 1, proportionnal, N);
V_Case1_Perf_Prop = Pricing_Phoenix(S0_values, Nmc, T, r, sigma, delta_t, Pi0, B_Ph_coeff_1, B_Y_coeff_1, B_Put_coeff_1, coeff_ph_1, coeff_y_1, 0, proportionnal, N);
V_Case2_Put_Prop = Pricing_Phoenix(S0_values, Nmc, T, r, sigma, delta_t, Pi0, B_Ph_coeff_2, B_Y_coeff_2, B_Put_coeff_2, coeff_ph_2, coeff_y_2, 1, proportionnal, N);
V_Case2_Perf_Prop = Pricing_Phoenix(S0_values, Nmc, T, r, sigma, delta_t, Pi0, B_Ph_coeff_2, B_Y_coeff_2, B_Put_coeff_2, coeff_ph_2, coeff_y_2, 0, proportionnal, N);

S0_index = find(S0_values == 100);
fprintf('Prices for S0 = 100 with Fixed Barriers:\n');
fprintf('Case 1 - Put Payoff: %.6f\n', V_Case1_Put_Prop(S0_index));
fprintf('Case 1 - Performance Payoff: %.6f\n', V_Case1_Perf_Prop(S0_index));
fprintf('Case 2 - Put Payoff: %.6f\n', V_Case2_Put_Prop(S0_index));
fprintf('Case 2 - Performance Payoff: %.6f\n', V_Case2_Perf_Prop(S0_index));
Graph_Pricing(S0_values, V_Case1_Put_Prop, V_Case1_Perf_Prop, V_Case2_Put_Prop, V_Case2_Perf_Prop, proportionnal);
% ---------- End proportional barriers 