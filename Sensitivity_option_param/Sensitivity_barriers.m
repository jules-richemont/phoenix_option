% Parameters
Pi0 = 1;            % Nominal value (100%)
S0 = 100;           % Fixed initial stock price
K = S0;             % Strike price
T = 5;              % Maturity in years
delta_t = 1;        % Observation frequency in years
sigma = 0.3;        % Volatility
r = 0.02;           % Risk-free rate
N = 10000;          % Number of simulations
timesteps = T / delta_t;

% Fixed coupons
C_Ph = 10;    % 10% of S0
C_Y = 5;      % 5% of S0

% Time vector
t = delta_t * (1:timesteps);

% Simulate N trajectories of S_t once for all analyses
Z = randn(N, timesteps);
S = zeros(N, timesteps);
S(:,1) = S0 * exp((r - 0.5 * sigma^2) * delta_t + sigma * sqrt(delta_t) .* Z(:,1));
for i = 2:timesteps
    S(:,i) = S(:,i-1) .* exp((r - 0.5 * sigma^2) * delta_t + sigma * sqrt(delta_t) .* Z(:,i));
end

%% Sensitivity to B_Ph
B_Ph_values = 100:2:140;    % Varying B_Ph from 100 to 140
B_Y_fixed = 80;             % Fixed B_Y
B_Put_fixed = 70;           % Fixed B_Put

V_BPh = zeros(length(B_Ph_values),1);

for idx = 1:length(B_Ph_values)
    B_Ph = B_Ph_values(idx);
    Payoff_Put = zeros(N,1);
    for n = 1:N
        PV = 0; autocall = false;
        for i = 1:timesteps
            if S(n,i) > B_Ph
                PV = (Pi0 + C_Ph) * exp(-r * t(i));
                autocall = true; break;
            elseif S(n,i) >= B_Y_fixed && S(n,i) <= B_Ph
                PV = PV + C_Y * exp(-r * t(i));
            end
        end
        if ~autocall
            S_N = S(n,end);
            if S_N > B_Ph
                PV = PV + (Pi0 + C_Ph) * exp(-r * T);
            elseif S_N >= B_Y_fixed && S_N <= B_Ph
                PV = PV + (Pi0 + C_Y) * exp(-r * T);
            elseif S_N > B_Put_fixed && S_N < B_Y_fixed
                PV = PV + Pi0 * exp(-r * T);
            else
                Payoff = max((K - S_N)/S0, 0);
                PV = PV + Payoff * exp(-r * T);
            end
        end
        Payoff_Put(n) = PV;
    end
    V_BPh(idx) = mean(Payoff_Put);
end

%% Sensitivity to B_Y
B_Y_values = 60:2:100;      % Varying B_Y from 60 to 100
B_Ph_fixed = 120;           % Fixed B_Ph
B_Put_fixed = 70;           % Fixed B_Put

V_BY = zeros(length(B_Y_values),1);

for idx = 1:length(B_Y_values)
    B_Y = B_Y_values(idx);
    Payoff_Put = zeros(N,1);
    for n = 1:N
        PV = 0; autocall = false;
        for i = 1:timesteps
            if S(n,i) > B_Ph_fixed
                PV = (Pi0 + C_Ph) * exp(-r * t(i));
                autocall = true; break;
            elseif S(n,i) >= B_Y && S(n,i) <= B_Ph_fixed
                PV = PV + C_Y * exp(-r * t(i));
            end
        end
        if ~autocall
            S_N = S(n,end);
            if S_N > B_Ph_fixed
                PV = PV + (Pi0 + C_Ph) * exp(-r * T);
            elseif S_N >= B_Y && S_N <= B_Ph_fixed
                PV = PV + (Pi0 + C_Y) * exp(-r * T);
            elseif S_N > B_Put_fixed && S_N < B_Y
                PV = PV + Pi0 * exp(-r * T);
            else
                Payoff = max((K - S_N)/S0, 0);
                PV = PV + Payoff * exp(-r * T);
            end
        end
        Payoff_Put(n) = PV;
    end
    V_BY(idx) = mean(Payoff_Put);
end

%% Sensitivity to B_Put
B_Put_values = 50:2:90;     % Varying B_Put from 50 to 90
B_Ph_fixed = 120;           % Fixed B_Ph
B_Y_fixed = 80;             % Fixed B_Y

V_BPut = zeros(length(B_Put_values),1);

for idx = 1:length(B_Put_values)
    B_Put = B_Put_values(idx);
    Payoff_Put = zeros(N,1);
    for n = 1:N
        PV = 0; autocall = false;
        for i = 1:timesteps
            if S(n,i) > B_Ph_fixed
                PV = (Pi0 + C_Ph) * exp(-r * t(i));
                autocall = true; break;
            elseif S(n,i) >= B_Y_fixed && S(n,i) <= B_Ph_fixed
                PV = PV + C_Y * exp(-r * t(i));
            end
        end
        if ~autocall
            S_N = S(n,end);
            if S_N > B_Ph_fixed
                PV = PV + (Pi0 + C_Ph) * exp(-r * T);
            elseif S_N >= B_Y_fixed && S_N <= B_Ph_fixed
                PV = PV + (Pi0 + C_Y) * exp(-r * T);
            elseif S_N > B_Put && S_N < B_Y_fixed
                PV = PV + Pi0 * exp(-r * T);
            else
                Payoff = max((K - S_N)/S0, 0);
                PV = PV + Payoff * exp(-r * T);
            end
        end
        Payoff_Put(n) = PV;
    end
    V_BPut(idx) = mean(Payoff_Put);
end

%% Plotting the results
figure;
subplot(1,3,1);
plot(B_Ph_values, V_BPh, 'b', 'LineWidth', 2);
xlabel('B_{Ph}');
ylabel('V(t=0, S_0, B_{Ph})');
title('Option Price vs B_{Ph}');
grid on;

subplot(1,3,2);
plot(B_Y_values, V_BY, 'r', 'LineWidth', 2);
xlabel('B_Y');
ylabel('V(t=0, S_0, B_Y)');
title('Option Price vs B_Y');
grid on;

subplot(1,3,3);
plot(B_Put_values, V_BPut, 'g', 'LineWidth', 2);
xlabel('B_{Put}');
ylabel('V(t=0, S_0, B_{Put})');
title('Option Price vs B_{Put}');
grid on;