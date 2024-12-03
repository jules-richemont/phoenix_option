% Parameters
Pi0 = 1;            % Nominal value (100%)
S0 = 100;           % Fixed initial stock price
K = S0;             % Strike price
delta_t = 1;        % Observation frequency in years
N = 100000;         % Number of simulations

% Common parameters
C_Ph = 0.1 * S0;    % Coupon at autocall (10% of S0)
C_Y = 0.05 * S0;    % Yearly coupon (5% of S0)

% Sensitivity ranges
sigma_values = 0.1:0.02:0.5;     % Volatility from 0.1 to 0.5
T_values = 1:1:10;               % Maturity from 1 to 10 years
r_values = 0:0.005:0.05;         % Risk-free rate from 0% to 5%

%% Sensitivity to Volatility sigma
V_sigma_caseA = zeros(length(sigma_values),1);
V_sigma_caseB = zeros(length(sigma_values),1);

for idx = 1:length(sigma_values)
    sigma = sigma_values(idx);
    T = 5;          % Fixed maturity
    r = 0.02;       % Fixed risk-free rate
    timesteps = T / delta_t;
    t = delta_t * (1:timesteps);
    
    % Barriers for Case a (proportional)
    B_Ph_a = 1.2 * S0;   % 120% of S0
    B_Y_a = 0.8 * S0;    % 80% of S0
    B_Put_a = 0.7 * S0;  % 70% of S0
    
    % Barriers for Case b (fixed)
    B_Ph_b = 120;        % Fixed at 120
    B_Y_b = 80;          % Fixed at 80
    B_Put_b = 70;        % Fixed at 70
    
    % Simulate asset paths
    Z = randn(N, timesteps);
    S = zeros(N, timesteps);
    S(:,1) = S0 * exp((r - 0.5 * sigma^2) * delta_t + sigma * sqrt(delta_t) .* Z(:,1));
    for i = 2:timesteps
        S(:,i) = S(:,i-1) .* exp((r - 0.5 * sigma^2) * delta_t + sigma * sqrt(delta_t) .* Z(:,i));
    end
    
    % Pricing for Case a
    Payoff_a = zeros(N,1);
    for n = 1:N
        PV = 0; autocall = false;
        for i = 1:timesteps
            if S(n,i) > B_Ph_a
                PV = (Pi0 + C_Ph) * exp(-r * t(i));
                autocall = true; break;
            elseif S(n,i) >= B_Y_a && S(n,i) <= B_Ph_a
                PV = PV + C_Y * exp(-r * t(i));
            end
        end
        if ~autocall
            S_N = S(n,end);
            if S_N > B_Ph_a
                PV = PV + (Pi0 + C_Ph) * exp(-r * T);
            elseif S_N >= B_Y_a && S_N <= B_Ph_a
                PV = PV + (Pi0 + C_Y) * exp(-r * T);
            elseif S_N > B_Put_a && S_N < B_Y_a
                PV = PV + Pi0 * exp(-r * T);
            else
                Payoff = max((K - S_N)/S0, 0);
                PV = PV + Payoff * exp(-r * T);
            end
        end
        Payoff_a(n) = PV;
    end
    V_sigma_caseA(idx) = mean(Payoff_a);
    
    % Pricing for Case b
    Payoff_b = zeros(N,1);
    for n = 1:N
        PV = 0; autocall = false;
        for i = 1:timesteps
            if S(n,i) > B_Ph_b
                PV = (Pi0 + C_Ph) * exp(-r * t(i));
                autocall = true; break;
            elseif S(n,i) >= B_Y_b && S(n,i) <= B_Ph_b
                PV = PV + C_Y * exp(-r * t(i));
            end
        end
        if ~autocall
            S_N = S(n,end);
            if S_N > B_Ph_b
                PV = PV + (Pi0 + C_Ph) * exp(-r * T);
            elseif S_N >= B_Y_b && S_N <= B_Ph_b
                PV = PV + (Pi0 + C_Y) * exp(-r * T);
            elseif S_N > B_Put_b && S_N < B_Y_b
                PV = PV + Pi0 * exp(-r * T);
            else
                Payoff = max((K - S_N)/S0, 0);
                PV = PV + Payoff * exp(-r * T);
            end
        end
        Payoff_b(n) = PV;
    end
    V_sigma_caseB(idx) = mean(Payoff_b);
end

%% Plotting V vs sigma
figure;
plot(sigma_values, V_sigma_caseA, 'b-', 'LineWidth', 2);
hold on;
plot(sigma_values, V_sigma_caseB, 'r--', 'LineWidth', 2);
xlabel('\sigma');
ylabel('V(t=0, S_0, \sigma)');
title('Option Price vs Volatility \sigma');
legend('Case a (Proportional Barriers)', 'Case b (Fixed Barriers)');
grid on;

%% Sensitivity to Maturity T
V_T_caseA = zeros(length(T_values),1);
V_T_caseB = zeros(length(T_values),1);

sigma = 0.3;      % Fixed volatility
r = 0.02;         % Fixed risk-free rate

for idx = 1:length(T_values)
    T = T_values(idx);
    timesteps = T / delta_t;
    t = delta_t * (1:timesteps);
    
    % Barriers for Case a
    B_Ph_a = 1.2 * S0;
    B_Y_a = 0.8 * S0;
    B_Put_a = 0.7 * S0;
    
    % Barriers for Case b
    B_Ph_b = 120;
    B_Y_b = 80;
    B_Put_b = 70;
    
    % Simulate asset paths
    Z = randn(N, timesteps);
    S = zeros(N, timesteps);
    S(:,1) = S0 * exp((r - 0.5 * sigma^2) * delta_t + sigma * sqrt(delta_t) .* Z(:,1));
    for i = 2:timesteps
        S(:,i) = S(:,i-1) .* exp((r - 0.5 * sigma^2) * delta_t + sigma * sqrt(delta_t) .* Z(:,i));
    end
    
    % Pricing for Case a
    Payoff_a = zeros(N,1);
    for n = 1:N
        PV = 0; autocall = false;
        for i = 1:timesteps
            t_i = t(i);
            if S(n,i) > B_Ph_a
                PV = (Pi0 + C_Ph) * exp(-r * t_i);
                autocall = true; break;
            elseif S(n,i) >= B_Y_a && S(n,i) <= B_Ph_a
                PV = PV + C_Y * exp(-r * t_i);
            end
        end
        if ~autocall
            S_N = S(n,end);
            if S_N > B_Ph_a
                PV = PV + (Pi0 + C_Ph) * exp(-r * T);
            elseif S_N >= B_Y_a && S_N <= B_Ph_a
                PV = PV + (Pi0 + C_Y) * exp(-r * T);
            elseif S_N > B_Put_a && S_N < B_Y_a
                PV = PV + Pi0 * exp(-r * T);
            else
                Payoff = max((K - S_N)/S0, 0);
                PV = PV + Payoff * exp(-r * T);
            end
        end
        Payoff_a(n) = PV;
    end
    V_T_caseA(idx) = mean(Payoff_a);
    
    % Pricing for Case b
    Payoff_b = zeros(N,1);
    for n = 1:N
        PV = 0; autocall = false;
        for i = 1:timesteps
            t_i = t(i);
            if S(n,i) > B_Ph_b
                PV = (Pi0 + C_Ph) * exp(-r * t_i);
                autocall = true; break;
            elseif S(n,i) >= B_Y_b && S(n,i) <= B_Ph_b
                PV = PV + C_Y * exp(-r * t_i);
            end
        end
        if ~autocall
            S_N = S(n,end);
            if S_N > B_Ph_b
                PV = PV + (Pi0 + C_Ph) * exp(-r * T);
            elseif S_N >= B_Y_b && S_N <= B_Ph_b
                PV = PV + (Pi0 + C_Y) * exp(-r * T);
            elseif S_N > B_Put_b && S_N < B_Y_b
                PV = PV + Pi0 * exp(-r * T);
            else
                Payoff = max((K - S_N)/S0, 0);
                PV = PV + Payoff * exp(-r * T);
            end
        end
        Payoff_b(n) = PV;
    end
    V_T_caseB(idx) = mean(Payoff_b);
end

%% Plotting V vs T
figure;
plot(T_values, V_T_caseA, 'b-', 'LineWidth', 2);
hold on;
plot(T_values, V_T_caseB, 'r--', 'LineWidth', 2);
xlabel('T (Maturity)');
ylabel('V(t=0, S_0, T)');
title('Option Price vs Maturity T');
legend('Case a (Proportional Barriers)', 'Case b (Fixed Barriers)');
grid on;

%% Sensitivity to Risk-free Rate r
V_r_caseA = zeros(length(r_values),1);
V_r_caseB = zeros(length(r_values),1);

sigma = 0.3;  % Fixed volatility
T = 5;        % Fixed maturity
timesteps = T / delta_t;
t = delta_t * (1:timesteps);

for idx = 1:length(r_values)
    r = r_values(idx);
    
    % Barriers for Case a
    B_Ph_a = 1.2 * S0;
    B_Y_a = 0.8 * S0;
    B_Put_a = 0.7 * S0;
    
    % Barriers for Case b
    B_Ph_b = 120;
    B_Y_b = 80;
    B_Put_b = 70;
    
    % Simulate asset paths
    Z = randn(N, timesteps);
    S = zeros(N, timesteps);
    S(:,1) = S0 * exp((r - 0.5 * sigma^2) * delta_t + sigma * sqrt(delta_t) .* Z(:,1));
    for i = 2:timesteps
        S(:,i) = S(:,i-1) .* exp((r - 0.5 * sigma^2) * delta_t + sigma * sqrt(delta_t) .* Z(:,i));
    end
    
    % Pricing for Case a
    Payoff_a = zeros(N,1);
    for n = 1:N
        PV = 0; autocall = false;
        for i = 1:timesteps
            t_i = t(i);
            if S(n,i) > B_Ph_a
                PV = (Pi0 + C_Ph) * exp(-r * t_i);
                autocall = true; break;
            elseif S(n,i) >= B_Y_a && S(n,i) <= B_Ph_a
                PV = PV + C_Y * exp(-r * t_i);
            end
        end
        if ~autocall
            S_N = S(n,end);
            if S_N > B_Ph_a
                PV = PV + (Pi0 + C_Ph) * exp(-r * T);
            elseif S_N >= B_Y_a && S_N <= B_Ph_a
                PV = PV + (Pi0 + C_Y) * exp(-r * T);
            elseif S_N > B_Put_a && S_N < B_Y_a
                PV = PV + Pi0 * exp(-r * T);
            else
                Payoff = max((K - S_N)/S0, 0);
                PV = PV + Payoff * exp(-r * T);
            end
        end
        Payoff_a(n) = PV;
    end
    V_r_caseA(idx) = mean(Payoff_a);
    
    % Pricing for Case b
    Payoff_b = zeros(N,1);
    for n = 1:N
        PV = 0; autocall = false;
        for i = 1:timesteps
            t_i = t(i);
            if S(n,i) > B_Ph_b
                PV = (Pi0 + C_Ph) * exp(-r * t_i);
                autocall = true; break;
            elseif S(n,i) >= B_Y_b && S(n,i) <= B_Ph_b
                PV = PV + C_Y * exp(-r * t_i);
            end
        end
        if ~autocall
            S_N = S(n,end);
            if S_N > B_Ph_b
                PV = PV + (Pi0 + C_Ph) * exp(-r * T);
            elseif S_N >= B_Y_b && S_N <= B_Ph_b
                PV = PV + (Pi0 + C_Y) * exp(-r * T);
            elseif S_N > B_Put_b && S_N < B_Y_b
                PV = PV + Pi0 * exp(-r * T);
            else
                Payoff = max((K - S_N)/S0, 0);
                PV = PV + Payoff * exp(-r * T);
            end
        end
        Payoff_b(n) = PV;
    end
    V_r_caseB(idx) = mean(Payoff_b);
end

%% Plotting V vs r
figure;
plot(r_values, V_r_caseA, 'b-', 'LineWidth', 2);
hold on;
plot(r_values, V_r_caseB, 'r--', 'LineWidth', 2);
xlabel('Risk-free Rate r');
ylabel('V(t=0, S_0, r)');
title('Option Price vs Risk-free Rate r');
legend('Case a (Proportional Barriers)', 'Case b (Fixed Barriers)');
grid on;