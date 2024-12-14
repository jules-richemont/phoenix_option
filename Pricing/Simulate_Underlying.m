function[S] = Simulate_Underlying(S0, r, sigma, T, N)
    % N: number of time steps
    % S: vector of simulated prices : size 1 x (N)
    S = zeros(1, N);
    S(1) = S0;
    dt = T/N;
    for i = 2:N
        S(i) = S(i-1) * exp((r - 0.5 * sigma^2) * (dt) + sigma * sqrt(dt) * randn);
    end
