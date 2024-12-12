function [V_Put, V_Perf] = Pricings_Phoenix(S0_values, N, T, r, sigma, delta_t, Pi0, B_Ph, B_Y, B_Put, coeff_ph, coeff_y)
    V_Put = zeros(length(S0_values), 1);
    V_Perf = zeros(length(S0_values), 1);

    for idx = 1:length(S0_values)
        S0 = S0_values(idx);
        K = S0;

        % Coupons
        C_Ph = coeff_ph * S0; % % of S0
        C_Y = coeff_y * S0; % % of S0


        % Initialize payoffs
        Payoff_Put = zeros(N, 1);
        Payoff_Perf = zeros(N, 1);

        timesteps = T / delta_t;
        % Version 1 
        S = Simulate_Trajectories(S0, r, sigma, T, timesteps, N);

        % Simulate trajectories and compute payoffs
        for n = 1:N
            % S = Simulate_Underlying(S0, r, sigma, T, T / delta_t); % Single trajectory

            % Compute Payoffs
            Payoff_Put(n) = Compute_Put_Payoff(S(n,:), K, Pi0, C_Ph, C_Y, T, r, B_Ph, B_Y, B_Put, delta_t);
            Payoff_Perf(n) = Compute_Perf_Payoff(S(n,:), Pi0, C_Ph, C_Y, T, r, B_Ph, B_Y, B_Put, delta_t);
        end

        % Average payoffs
        V_Put(idx) = mean(Payoff_Put);
        V_Perf(idx) = mean(Payoff_Perf);
    end
end
