% Parameters
Pi0 = 1;            
S0 = 100;           
K = S0;             
T = 5;              
delta_t = 1;        
sigma = 0.3;       

N = 10000;         
timesteps = T / delta_t;

% Fixed barriers
B_Ph = 120;    
B_Y = 80;      
B_Put = 70;    

% Time vector
t = delta_t * (1:timesteps);

% Simulate N trajectories of S_t once for both analyses
Z = randn(N, timesteps);
S = zeros(N, timesteps);
S(:,1) = S0;
for i = 2:timesteps
    S(:,i) = S(:,i-1) .* exp((r - 0.5 * sigma^2) * delta_t + sigma * sqrt(delta_t) .* Z(:,i));
end

%% Sensitivity to C_Ph
C_Ph_values = 0:0.5:20;    
C_Y_fixed = 5;             

V_CPh = zeros(length(C_Ph_values),1);

for idx = 1:length(C_Ph_values)
    C_Ph = C_Ph_values(idx);
    Payoff_Put = zeros(N,1);
    for n = 1:N
        PV = 0; autocall = false;
        for i = 1:timesteps
            if S(n,i) > B_Ph
                PV = (Pi0 + C_Ph) * exp(-r * t(i));
                autocall = true; break;
            elseif S(n,i) >= B_Y && S(n,i) <= B_Ph
                PV = PV + C_Y_fixed * exp(-r * t(i));
            end
        end
        if ~autocall
            S_N = S(n,end);
            if S_N > B_Ph
                PV = PV + (Pi0 + C_Ph) * exp(-r * T);
            elseif S_N >= B_Y && S_N <= B_Ph
                PV = PV + (Pi0 + C_Y_fixed) * exp(-r * T);
            elseif S_N > B_Put && S_N < B_Y
                PV = PV + Pi0 * exp(-r * T);
            else
                Payoff = max((K - S_N)/S0, 0);
                PV = PV + Payoff * exp(-r * T);
            end
        end
        Payoff_Put(n) = PV;
    end
    V_CPh(idx) = mean(Payoff_Put);
end

%% Sensitivity to C_Y
C_Y_values = 0:0.5:10;     
C_Ph_fixed = 10;           

V_CY = zeros(length(C_Y_values),1);

for idx = 1:length(C_Y_values)
    C_Y = C_Y_values(idx);
    Payoff_Put = zeros(N,1);
    for n = 1:N
        PV = 0; autocall = false;
        for i = 1:timesteps
            if S(n,i) > B_Ph
                PV = (Pi0 + C_Ph_fixed) * exp(-r * t(i));
                autocall = true; break;
            elseif S(n,i) >= B_Y && S(n,i) <= B_Ph
                PV = PV + C_Y * exp(-r * t(i));
            end
        end
        if ~autocall
            S_N = S(n,end);
            if S_N > B_Ph
                PV = PV + (Pi0 + C_Ph_fixed) * exp(-r * T);
            elseif S_N >= B_Y && S_N <= B_Ph
                PV = PV + (Pi0 + C_Y) * exp(-r * T);
            elseif S_N > B_Put && S_N < B_Y
                PV = PV + Pi0 * exp(-r * T);
            else
                Payoff = max((K - S_N)/S0, 0);
                PV = PV + Payoff * exp(-r * T);
            end
        end
        Payoff_Put(n) = PV;
    end
    V_CY(idx) = mean(Payoff_Put);
end

%% Plotting the results
figure;
subplot(1,2,1);
plot(C_Ph_values, V_CPh, 'b', 'LineWidth', 2);
xlabel('C_{Ph}');
ylabel('V(t=0, S_0, C_{Ph})');
title('Option Price vs C_{Ph} (C_Y = 5)');
grid on;

subplot(1,2,2);
plot(C_Y_values, V_CY, 'r', 'LineWidth', 2);
xlabel('C_Y');
ylabel('V(t=0, S_0, C_Y)');
title('Option Price vs C_Y (C_{Ph} = 10)');
grid on;