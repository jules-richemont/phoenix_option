% Parameters
Pi0 = 1;            % Nominal value (100%)
T = 10;             % Maturity in years
delta_t = 1;        % Observation frequency in years
sigma = 0.3;        % Volatility
r = 0.02;           % Risk-free rate
N = 100000;         % Number of simulations
timesteps = T / delta_t;

S0_values = 80:1:120;   % S0 varying from 80 to 120

% Fixed barriers
B_Ph = 120;    % Fixed autocall barrier
B_Y = 80;      % Fixed coupon barrier
B_Put = 70;    % Fixed put barrier

% Pre-allocate arrays to store prices
V_Case1_Put = zeros(length(S0_values),1);
V_Case1_Perf = zeros(length(S0_values),1);
V_Case2_Put = zeros(length(S0_values),1);
V_Case2_Perf = zeros(length(S0_values),1);

% Loop over S0 values
for idx = 1:length(S0_values)
    S0 = S0_values(idx);
    K = S0;
    
    % Coupons for Case 1
    C_Ph1 = 0.1 * S0;   % 10% of S0
    C_Y1 = 0.05 * S0;   % 5% of S0
    
    % Coupons for Case 2
    C_Ph2 = 0.08 * S0;  % 8% of S0
    C_Y2 = 0.05 * S0;   % 5% of S0
    
    % Time vector
    t = delta_t * (1:timesteps);
    
    % Simulate N trajectories of S_t
    Z = randn(N, timesteps); % N x timesteps matrix of standard normal random variables
    S = zeros(N, timesteps);
    S(:,1) = S0;
    for i = 2:timesteps
        S(:,i) = S(:,i-1) .* exp((r - 0.5 * sigma^2) * delta_t + sigma * sqrt(delta_t) * Z(:,i));
    end
    
    %% Pricing for Case 1 with Put Payoff
    Payoff_Case1_Put = zeros(N,1);
    for n = 1:N
        PV = 0; % Present Value for this trajectory
        autocall = false;
        for i = 1:timesteps
            if S(n,i) > B_Ph
                % Autocall triggered
                PV = (Pi0 + C_Ph1) * exp(-r * t(i));
                autocall = true;
                break;
            elseif S(n,i) >= B_Y && S(n,i) <= B_Ph
                % Coupon payment
                PV = PV + C_Y1 * exp(-r * t(i));
            end
        end
        if ~autocall
            % Evaluate payoff at maturity
            S_N = S(n,end);
            if S_N > B_Ph
                PV = PV + (Pi0 + C_Ph1) * exp(-r * T);
            elseif S_N >= B_Y && S_N <= B_Ph
                PV = PV + (Pi0 + C_Y1) * exp(-r * T);
            elseif S_N > B_Put && S_N < B_Y
                PV = PV + Pi0 * exp(-r * T);
            elseif S_N <= B_Put
                % Put Payoff
                Payoff = max((K - S_N)/S0, 0);
                PV = PV + Payoff * exp(-r * T);
            end
        end
        Payoff_Case1_Put(n) = PV;
    end
    V_Case1_Put(idx) = mean(Payoff_Case1_Put);
    
    %% Pricing for Case 1 with Performance Payoff
    Payoff_Case1_Perf = zeros(N,1);
    for n = 1:N
        PV = 0;
        autocall = false;
        for i = 1:timesteps
            if S(n,i) > B_Ph
                PV = (Pi0 + C_Ph1) * exp(-r * t(i));
                autocall = true;
                break;
            elseif S(n,i) >= B_Y && S(n,i) <= B_Ph
                PV = PV + C_Y1 * exp(-r * t(i));
            end
        end
        if ~autocall
            S_N = S(n,end);
            if S_N > B_Ph
                PV = PV + (Pi0 + C_Ph1) * exp(-r * T);
            elseif S_N >= B_Y && S_N <= B_Ph
                PV = PV + (Pi0 + C_Y1) * exp(-r * T);
            elseif S_N > B_Put && S_N < B_Y
                PV = PV + Pi0 * exp(-r * T);
            elseif S_N <= B_Put
                % Performance Payoff
                Payoff = S_N / S0;
                PV = PV + Payoff * exp(-r * T);
            end
        end
        Payoff_Case1_Perf(n) = PV;
    end
    V_Case1_Perf(idx) = mean(Payoff_Case1_Perf);
    
    %% Pricing for Case 2 with Put Payoff
    Payoff_Case2_Put = zeros(N,1);
    for n = 1:N
        PV = 0; % Present Value for this trajectory
        autocall = false;
        for i = 1:timesteps
            if S(n,i) > B_Ph
                % Autocall triggered
                PV = (Pi0 + C_Ph2) * exp(-r * t(i));
                autocall = true;
                break;
            elseif S(n,i) >= B_Y && S(n,i) <= B_Ph
                % Coupon payment
                PV = PV + C_Y2 * exp(-r * t(i));
            end
        end
        if ~autocall
            % Evaluate payoff at maturity
            S_N = S(n,end);
            if S_N > B_Ph
                PV = PV + (Pi0 + C_Ph2) * exp(-r * T);
            elseif S_N >= B_Y && S_N <= B_Ph
                PV = PV + (Pi0 + C_Y2) * exp(-r * T);
            elseif S_N > B_Put && S_N < B_Y
                PV = PV + Pi0 * exp(-r * T);
            elseif S_N <= B_Put
                % Put Payoff
                Payoff = max((K - S_N)/S0, 0);
                PV = PV + Payoff * exp(-r * T);
            end
        end
        Payoff_Case2_Put(n) = PV;
    end
    V_Case2_Put(idx) = mean(Payoff_Case2_Put);
    
    %% Pricing for Case 2 with Performance Payoff
    Payoff_Case2_Perf = zeros(N,1);
    for n = 1:N
        PV = 0; % Present Value for this trajectory
        autocall = false;
        for i = 1:timesteps
            if S(n,i) > B_Ph
                % Autocall triggered
                PV = (Pi0 + C_Ph2) * exp(-r * t(i));
                autocall = true;
                break;
            elseif S(n,i) >= B_Y && S(n,i) <= B_Ph
                % Coupon payment
                PV = PV + C_Y2 * exp(-r * t(i));
            end
        end
        if ~autocall
            % Evaluate payoff at maturity
            S_N = S(n,end);
            if S_N > B_Ph
                PV = PV + (Pi0 + C_Ph2) * exp(-r * T);
            elseif S_N >= B_Y && S_N <= B_Ph
                PV = PV + (Pi0 + C_Y2) * exp(-r * T);
            elseif S_N > B_Put && S_N < B_Y
                PV = PV + Pi0 * exp(-r * T);
            elseif S_N <= B_Put
                % Performance Payoff
                Payoff = S_N / S0;
                PV = PV + Payoff * exp(-r * T);
            end
        end
        Payoff_Case2_Perf(n) = PV;
    end
    V_Case2_Perf(idx) = mean(Payoff_Case2_Perf);
    
end

% Plotting the results
figure;
plot(S0_values, V_Case1_Put, 'b', 'LineWidth', 2);
hold on;
plot(S0_values, V_Case1_Perf, 'b--', 'LineWidth', 2);
plot(S0_values, V_Case2_Put, 'r', 'LineWidth', 2);
plot(S0_values, V_Case2_Perf, 'r--', 'LineWidth', 2);
xlabel('S_0');
ylabel('V(S_0, t_0)');
title('Option Prices vs Initial Stock Price S_0 with Fixed Barriers');
legend('Case 1 - Put Payoff', 'Case 1 - Performance Payoff', 'Case 2 - Put Payoff', 'Case 2 - Performance Payoff');
grid on;

% Display prices for S0 = 100
S0_index = find(S0_values == 100);
fprintf('Prices for S0 = 100 with Fixed Barriers:\n');
fprintf('Case 1 - Put Payoff: %f\n', V_Case1_Put(S0_index));
fprintf('Case 1 - Performance Payoff: %f\n', V_Case1_Perf(S0_index));
fprintf('Case 2 - Put Payoff: %f\n', V_Case2_Put(S0_index));
fprintf('Case 2 - Performance Payoff: %f\n', V_Case2_Perf(S0_index));