function Graph_Delta(Nmc, K, r, sigma, B, T, S_min, S_max, h)
    S0_values = linspace(S_min, S_max, 200);
    Delta_values = zeros(size(S0_values));
    
    for i = 1:length(S0_values)
        S0 = S0_values(i);
        Delta_values(i) = Compute_Delta(Nmc, S0, K, r, sigma, B, T, h);
    end
    
    % Tracer le graphe
    figure;
    plot(S0_values, Delta_values, '-o');
    title('Delta vs S_0 for Down-and-In Put Option');
    xlabel('S_0');
    ylabel('\Delta(S_0)');
    grid on;
end
