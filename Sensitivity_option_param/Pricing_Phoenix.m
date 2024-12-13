function [V] = Pricing_Phoenix(S0_values, N, T, r, sigma, delta_t, Pi0, B_Ph, B_Y, B_Put, coeff_ph, coeff_y, put_or_perf)
    % if put_or_perf == 1 : put
    % if put_or_perf == 0 : perf
    V = zeros(length(S0_values), 1);


    for idx = 1:length(S0_values)
        S0 = S0_values(idx);
        K = S0;

        % Coupons
        C_Ph = coeff_ph * S0; % % of S0
        C_Y = coeff_y * S0; % % of S0


        % Initialize payoffs
        Payoff = zeros(N, 1);


        timesteps = T / delta_t;
        % Version 1 
        S = Simulate_Trajectories(S0, r, sigma, T, timesteps, N);

        % Simulate trajectories and compute payoffs
        for n = 1:N
            % S = Simulate_Underlying(S0, r, sigma, T, T / delta_t); % Single trajectory

            % Compute Payoffs
            Payoff(n) = Compute_Payoff(S(n,:), K, Pi0, C_Ph, C_Y, T, r, B_Ph, B_Y, B_Put, delta_t, put_or_perf);
        end

        % Average payoff
        V(idx) = mean(Payoff);
    end
end
