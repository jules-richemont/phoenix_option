function[price] = Pricing_Barriere_Put_Down_In(Nmc, S0, K, r, sigma, B, T, N)
    payoff = zeros(1, Nmc);
    for n = 1:Nmc
        S = Simulate_Underlying(S0, r, sigma, T, N);
        % Check if the barrier is crossed : if the option is activated
        if min(S) < B
            % Compute the payoff
            payoff(n) = max(K - S(end), 0);
        else
            payoff(n) = 0;
        end
    end
    price = mean(payoff)*exp(-r*T);
end