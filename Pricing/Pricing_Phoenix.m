function [V] = Pricing_Phoenix(S0_values, Nmc, T, r, sigma, delta_t_observation, Pi0, B_Ph, B_Y, B_Put, coeff_ph, coeff_y, put_or_perf, proportionnal, N)
    % if put_or_perf == 1 : put
    % if put_or_perf == 0 : perf
    % be carefull : if prportionnal == 1, the barriers are proportionnal to the initial value S0 : B_Ph, B_Y, B_Put are coefficients
    % if prportionnal == 0, the barriers are fixed values and B_Ph, B_Y, B_Put are the values of the barriers
    V = zeros(length(S0_values), 1);
    % Coupons
    C_Ph = coeff_ph * Pi0; % % of S0
    C_Y = coeff_y * Pi0; % % of S0
    % Store original barrier values
    original_B_Ph = B_Ph;
    original_B_Y = B_Y;
    original_B_Put = B_Put;
    for idx = 1:length(S0_values)
        S0 = S0_values(idx);
        K = S0;

        % Initialize payoffs
        Payoff = zeros(Nmc, 1);
        % Version 1 
        S = Simulate_Trajectories(S0, r, sigma, T, N, Nmc);
        % print la taille d'une ligne de S



        % Simulate trajectories and compute payoffs
        for n = 1:Nmc
            % S = Simulate_Underlying(S0, r, sigma, T, T / delta_t); % Single trajectory
            if proportionnal == 1
                % fprintf('S0: %.2f, B_Ph: %.2f, B_Y: %.2f, B_Put: %.2f\n', S0, B_Ph, B_Y, B_Put);

                % It means that the barrier is proportionnal to the initial value S0
                B_Ph = original_B_Ph * S(n, 1);
                B_Y = original_B_Y * S(n, 1);
                B_Put = original_B_Put * S(n, 1);
            end
            % Compute Payoffs
            Payoff(n) = Compute_Payoff(S(n,:), K, Pi0, C_Ph, C_Y, T, r, B_Ph, B_Y, B_Put, delta_t_observation, put_or_perf);
            % if S(n, 1) == 100
            %     fprintf('S0: %.2f, B_Ph: %.2f, B_Y: %.2f, B_Put: %.2f\n', S0, B_Ph, B_Y, B_Put);
            %     fprintf('Payoff: %.2f\n', Payoff(n));
            % end
            
        end

        % Average payoff
        V(idx) = mean(Payoff);
    end
end
