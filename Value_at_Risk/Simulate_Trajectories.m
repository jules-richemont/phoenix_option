function S = Simulate_Trajectories(S0, r, sigma, T, timesteps, N)
    dt = T / timesteps;
    Z = randn(N, timesteps); % Générer toutes les variables aléatoires en une fois
    S = zeros(N, timesteps);
    S(:,1) = S0; % Initialisation
    for i = 2:timesteps
        S(:,i) = S(:,i-1) .* exp((r - 0.5 * sigma^2) * dt + sigma * sqrt(dt) * Z(:,i));
    end
end
