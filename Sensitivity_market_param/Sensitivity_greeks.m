% Parameters
Pi0 = 1;    % Nominal value (100%)
T = 5; % Maturity in years
delta_t = 1; % Observation frequency in years
r = 0.02;   % Risk-free rate
N = 100000; % Number of simulations
timesteps = T / delta_t;

% Fixed barriers
B_Ph = 120;% Fixed autocall barrier
B_Y = 80;% Fixed coupon barrier
B_Put = 70;% Fixed put barrier

% Coupons
C_Ph = 0.1 * 100;% 10% of S0 = 100
C_Y = 0.05 * 100; % 5% of S0 = 100

% Small increment h
h = 0.1;

% S0 values for plotting Delta and Gamma
S0_values = 80:1:120;

% Pre-allocate arrays for Greeks
Delta = zeros(length(S0_values),1);
Gamma = zeros(length(S0_values),1);
Vega_S0 = zeros(length(S0_values),1);

% Fixed volatility for Delta and Gamma
sigma = 0.3;

% Loop over S0 values
for idx = 1:length(S0_values)
    S0 = S0_values(idx);
    K = S0;  % Strike price
    
    % Simulate common random numbers
    Z = randn(N, timesteps);
    
    %% Compute V(S0 + h)
    S_plus = zeros(N, timesteps);
    S_plus(:,1) = (S0 + h) * exp((r - 0.5 * sigma^2) * delta_t + sigma * sqrt(delta_t) .* Z(:,1));
    for i = 2:timesteps
        S_plus(:,i) = S_plus(:,i-1) .* exp((r - 0.5 * sigma^2) * delta_t + sigma * sqrt(delta_t) .* Z(:,i));
    end
    V_plus = optionPricing(S_plus, B_Ph, B_Y, B_Put, C_Ph, C_Y, K, r, T, delta_t, Pi0);
    
    %% Compute V(S0 - h)
    S_minus = zeros(N, timesteps);
    S_minus(:,1) = (S0 - h) * exp((r - 0.5 * sigma^2) * delta_t + sigma * sqrt(delta_t) .* Z(:,1));
    for i = 2:timesteps
        S_minus(:,i) = S_minus(:,i-1) .* exp((r - 0.5 * sigma^2) * delta_t + sigma * sqrt(delta_t) .* Z(:,i));
    end
    V_minus = optionPricing(S_minus, B_Ph, B_Y, B_Put, C_Ph, C_Y, K, r, T, delta_t, Pi0);
    
    %% Compute V(S0)
    S = zeros(N, timesteps);
    S(:,1) = S0 * exp((r - 0.5 * sigma^2) * delta_t + sigma * sqrt(delta_t) .* Z(:,1));
    for i = 2:timesteps
        S(:,i) = S(:,i-1) .* exp((r - 0.5 * sigma^2) * delta_t + sigma * sqrt(delta_t) .* Z(:,i));
    end
    V = optionPricing(S, B_Ph, B_Y, B_Put, C_Ph, C_Y, K, r, T, delta_t, Pi0);
    
    %% Calculate Delta and Gamma
    Delta(idx) = (mean(V_plus) - mean(V_minus)) / (2 * h);
    Gamma(idx) = (mean(V_plus) - 2 * mean(V) + mean(V_minus)) / (h^2);
end

%% Plot Delta and Gamma vs S0
figure;
subplot(1,2,1);
plot(S0_values, Delta, 'b', 'LineWidth', 2);
xlabel('S_0');
ylabel('\Delta');
title('Delta vs S_0');
grid on;

subplot(1,2,2);
plot(S0_values, Gamma, 'r', 'LineWidth', 2);
xlabel('S_0');
ylabel('\Gamma');
title('Gamma vs S_0');
grid on;

%% Vega calculation for varying S0
% Small increment h for sigma
h_sigma = 0.01;
sigma_base = 0.3;
Vega_S0 = zeros(length(S0_values),1);

% Simulate common random numbers
Z = randn(N, timesteps);

for idx = 1:length(S0_values)
    S0 = S0_values(idx);
    K = S0;
    
    % Compute V(sigma + h)
    sigma_plus = sigma_base + h_sigma;
    S_sigma_plus = simulateAssetPaths(S0, sigma_plus, r, delta_t, N, timesteps, Z);
    V_sigma_plus = optionPricing(S_sigma_plus, B_Ph, B_Y, B_Put, C_Ph, C_Y, K, r, T, delta_t, Pi0);
    
    % Compute V(sigma - h)
    sigma_minus = sigma_base - h_sigma;
    S_sigma_minus = simulateAssetPaths(S0, sigma_minus, r, delta_t, N, timesteps, Z);
    V_sigma_minus = optionPricing(S_sigma_minus, B_Ph, B_Y, B_Put, C_Ph, C_Y, K, r, T, delta_t, Pi0);
    
    % Vega calculation
    Vega_S0(idx) = (mean(V_sigma_plus) - mean(V_sigma_minus)) / (2 * h_sigma);
end

%% Plot Vega vs S0
figure;
plot(S0_values, Vega_S0, 'g', 'LineWidth', 2);
xlabel('S_0');
ylabel('Vega');
title('Vega vs S_0');
grid on;

%% Vega calculation for varying sigma at S0 = 100
sigma_values = 0.1:0.02:0.5;
Vega_sigma = zeros(length(sigma_values),1);

S0 = 100;
K = S0;

% Simulate common random numbers
Z = randn(N, timesteps);

for idx = 1:length(sigma_values)
    sigma_base = sigma_values(idx);
    
    % Compute V(sigma + h)
    sigma_plus = sigma_base + h_sigma;
    S_sigma_plus = simulateAssetPaths(S0, sigma_plus, r, delta_t, N, timesteps, Z);
    V_sigma_plus = optionPricing(S_sigma_plus, B_Ph, B_Y, B_Put, C_Ph, C_Y, K, r, T, delta_t, Pi0);
    
    % Compute V(sigma - h)
    sigma_minus = sigma_base - h_sigma;
    S_sigma_minus = simulateAssetPaths(S0, sigma_minus, r, delta_t, N, timesteps, Z);
    V_sigma_minus = optionPricing(S_sigma_minus, B_Ph, B_Y, B_Put, C_Ph, C_Y, K, r, T, delta_t, Pi0);
    
    % Vega calculation
    Vega_sigma(idx) = (mean(V_sigma_plus) - mean(V_sigma_minus)) / (2 * h_sigma);
end

%% Plot Vega vs sigma
figure;
plot(sigma_values, Vega_sigma, 'm', 'LineWidth', 2);
xlabel('\sigma');
ylabel('Vega');
title('Vega vs \sigma at S_0 = 100');
grid on;

%% Function to simulate asset paths
function S = simulateAssetPaths(S0, sigma, r, delta_t, N, timesteps, Z)
    S = zeros(N, timesteps);
    S(:,1) = S0 * exp((r - 0.5 * sigma^2) * delta_t + sigma * sqrt(delta_t) .* Z(:,1));
    for i = 2:timesteps
        S(:,i) = S(:,i-1) .* exp((r - 0.5 * sigma^2) * delta_t + sigma * sqrt(delta_t) .* Z(:,i));
    end
end

%% Function to calculate option price
function V = optionPricing(S, B_Ph, B_Y, B_Put, C_Ph, C_Y, K, r, T, delta_t, Pi0)
    N = size(S,1);
    timesteps = size(S,2);
    t = delta_t * (1:timesteps);
    V = zeros(N,1);
    for n = 1:N
        PV = 0; autocall = false;
        for i = 1:timesteps
            if S(n,i) > B_Ph
                PV = (Pi0 + C_Ph) * exp(-r * t(i));
                autocall = true; break;
            elseif S(n,i) >= B_Y && S(n,i) <= B_Ph
                PV = PV + C_Y * exp(-r * t(i));
            end
        end
        if ~autocall
            S_N = S(n,end);
            if S_N > B_Ph
                PV = PV + (Pi0 + C_Ph) * exp(-r * T);
            elseif S_N >= B_Y && S_N <= B_Ph
                PV = PV + (Pi0 + C_Y) * exp(-r * T);
            elseif S_N > B_Put && S_N < B_Y
                PV = PV + Pi0 * exp(-r * T);
            else
                % Put Payoff
                Payoff = max((K - S_N)/K, 0);
                PV = PV + Payoff * exp(-r * T);
            end
        end
        V(n) = PV;
    end
end