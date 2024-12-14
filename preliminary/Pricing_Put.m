function[price]  = Pricing_Put(Nmc, S0, K, r, sigma, T) 
    payoff = zeros(1, Nmc);
    for n = 1:Nmc
        S = Simulate_Underlying(S0, r, sigma, T, 100);
        payoff(n) = max(K - S(end), 0);
    end
    price = mean(payoff)*exp(-r*T);
end