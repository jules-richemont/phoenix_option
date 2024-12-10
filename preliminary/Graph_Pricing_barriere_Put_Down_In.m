function[] = Graph_Pricing_barriere_Put_Down_In(Nmc, K, r, sigma, B, T, S_min, S_max)
    
    % COmpute the graph S_0 -> V(t=0, S_0)
    S0 = linspace(S_min, S_max, 100);
    price = zeros(1, 100);
    for i = 1:100
        price(i) = Pricing_Barriere_Put_Down_In(Nmc, S0(i), K, r, sigma, B, T);
    end

    % Plot the function
    plot(S0, price, 'r');
    title('Price of a Down-and-In Barrier Put Option');
    xlabel('S_0');
    ylabel('Price');
    % Hold off to finish the plot
    hold off;

end


