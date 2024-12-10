function [Delta] = Compute_Delta(Nmc, S0, K, r, sigma, B, T, h)
    % - h: small perturbation for finite differences
    
    payoff_plus = zeros(1, Nmc);
    payoff_minus = zeros(1, Nmc);
    N = 100;
    for n = 1:Nmc
        g = randn(1, N); % shared random numbers
        
        S_plus = Simulate_Underlying_Shared(S0 + h, r, sigma, T, N, g);
        S_minus = Simulate_Underlying_Shared(S0 - h, r, sigma, T, N, g);
        
    
        if min(S_plus) < B
            payoff_plus(n) = max(K - S_plus(end), 0);
        else
            payoff_plus(n) = 0;
        end
        if min(S_minus) < B
            payoff_minus(n) = max(K - S_minus(end), 0);
        else
            payoff_minus(n) = 0;
        end
    end
    
    %We use finite difference
    Delta = exp(-r * T) * (mean(payoff_plus) - mean(payoff_minus)) / (2 * h);
end
