% Parameters
Pi0 = 1;            % Nominal value (100%)
S0 = 100;           % Initial stock price
K = S0;             % Strike price
T = 5;              % Maturity in years
delta_t = 1;        % Observation frequency in years
sigma = 0.3;        % Volatility
r = 0.02;           % Risk-free rate
N = 100000;         % Number of simulations
timesteps = T / delta_t;

% Fixed barriers
B_Ph = 120;    % Fixed autocall barrier
B_Y = 80;      % Fixed coupon barrier
B_Put = 70;    % Fixed put barrier

% Coupons
C_Ph = 0.1 * S0;   % 10% of S0
C_Y = 0.05 * S0;   % 5% of S0

% Simulate asset paths
% Simulate asset paths
S = Simulate_Trajectories(S0, r, sigma, delta_t, timesteps, N);

% Time vector
t = delta_t * (1:timesteps);

% Simulate payoffs without discounting for Put Payoff
Payoff_Put = zeros(N,1);
for n = 1:N
    autocall = false;
    for i = 1:timesteps
        if S(n,i) > B_Ph
            % Autocall triggered
            Payoff_Put(n) = Pi0 + C_Ph;
            autocall = true;
            break;
        end
    end
    if ~autocall
        S_N = S(n,end);
        if S_N > B_Ph
            Payoff_Put(n) = Pi0 + C_Ph;
        elseif S_N >= B_Y && S_N <= B_Ph
            Payoff_Put(n) = Pi0 + C_Y;
        elseif S_N > B_Put && S_N < B_Y
            Payoff_Put(n) = Pi0;
        else
            % Put Payoff
            Payoff_Put(n) = max((K - S_N)/S0, 0);
        end
    end
end

% Simulate payoffs without discounting for Performance Payoff
Payoff_Perf = zeros(N,1);
for n = 1:N
    autocall = false;
    for i = 1:timesteps
        if S(n,i) > B_Ph
            % Autocall triggered
            Payoff_Perf(n) = Pi0 + C_Ph;
            autocall = true;
            break;
        end
    end
    if ~autocall
        S_N = S(n,end);
        if S_N > B_Ph
            Payoff_Perf(n) = Pi0 + C_Ph;
        elseif S_N >= B_Y && S_N <= B_Ph
            Payoff_Perf(n) = Pi0 + C_Y;
        elseif S_N > B_Put && S_N < B_Y
            Payoff_Perf(n) = Pi0;
        else
            % Performance Payoff
            Payoff_Perf(n) = S_N / S0;
        end
    end
end

% Calculate V0 (option price at t = 0) for both payoffs
V0_Put = mean(Payoff_Put .* exp(-r * T));
V0_Perf = mean(Payoff_Perf .* exp(-r * T));

% For VaR calculation, we need V_T - V0 (without discounting)
Losses_Put = Payoff_Put - V0_Put;
Losses_Perf = Payoff_Perf - V0_Perf;

% Sort the losses for empirical CDF
[Losses_Put_sorted, idx_Put] = sort(Losses_Put);
F_Put = (1:N)' / N;

[Losses_Perf_sorted, idx_Perf] = sort(Losses_Perf);
F_Perf = (1:N)' / N;

% Plot empirical CDFs
figure;
subplot(1,2,1);
plot(Losses_Put_sorted, F_Put, 'b', 'LineWidth', 2);
xlabel('V_T - V_0');
ylabel('F(V_T - V_0)');
title('CDF of V_T - V_0 for Put Payoff');
grid on;

subplot(1,2,2);
plot(Losses_Perf_sorted, F_Perf, 'r', 'LineWidth', 2);
xlabel('V_T - V_0');
ylabel('F(V_T - V_0)');
title('CDF of V_T - V_0 for Performance Payoff');
grid on;

% Calculate VaR at alpha = 10% and alpha = 1% using ordering algorithm
alphas = [0.10, 0.01];

% VaR for Put Payoff
VaR_Put = zeros(length(alphas),1);
for idx = 1:length(alphas)
    VaR_Put(idx) = Quantile(Losses_Put, alphas(idx));
end

% VaR for Performance Payoff
VaR_Perf = zeros(length(alphas),1);
for idx = 1:length(alphas)
    VaR_Perf(idx) = Quantile(Losses_Perf, alphas(idx));
end

% Mark VaR on CDF plots
subplot(1,2,1);
hold on;
for idx = 1:length(alphas)
    % Find the index closest to alpha in the empirical CDF
    [~, alpha_idx] = min(abs(F_Put - alphas(idx)));
    plot(Losses_Put_sorted(alpha_idx), F_Put(alpha_idx), 'ko', 'MarkerFaceColor', 'k');
    text(Losses_Put_sorted(alpha_idx), F_Put(alpha_idx), sprintf(' VaR at %.0f%%: %.4f', alphas(idx)*100, Losses_Put_sorted(alpha_idx)), 'VerticalAlignment', 'bottom');
end

subplot(1,2,2);
hold on;
for idx = 1:length(alphas)
    % Find the index closest to alpha in the empirical CDF
    [~, alpha_idx] = min(abs(F_Perf - alphas(idx)));
    plot(Losses_Perf_sorted(alpha_idx), F_Perf(alpha_idx), 'ko', 'MarkerFaceColor', 'k');
    text(Losses_Perf_sorted(alpha_idx), F_Perf(alpha_idx), sprintf(' VaR at %.0f%%: %.4f', alphas(idx)*100, Losses_Perf_sorted(alpha_idx)), 'VerticalAlignment', 'bottom');
end
