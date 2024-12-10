function [S] = Simulate_Underlying_Shared(S0, r, sigma, T, N, g)
    % simulate underlying with predefined random ran
    % - g: shared random numbers (1 x N)
    % Output:
    % - S: vector of simulated prices (1 x (N+1))
    S = zeros(1, N+1);
    S(1) = S0;
    dt = T / N;
    for i = 1:N
        S(i+1) = S(i) * exp((r - 0.5 * sigma^2) * dt + sigma * sqrt(dt) * g(i));
    end
end
