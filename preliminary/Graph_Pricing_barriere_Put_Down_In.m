function[] = Graph_Pricing_barriere_Put_Down_In(Nmc, K, r, sigma, B, T, S_min, S_max, N)
    
    % COmpute the graph 
    S0 = linspace(S_min, S_max, 100);
    price = zeros(1, 100);
    price_put = zeros(1, 100);
    for i = 1:100
        price(i) = Pricing_Barriere_Put_Down_In(Nmc, S0(i), K, r, sigma, B, T, N);
        price_put(i) = Pricing_Put(Nmc, S0(i), K, r, sigma, T);
    end

    % Plot
    figure;
    plot(S0, price, 'r', 'LineWidth', 2);
    hold on;
    plot(S0, price_put, 'b', 'LineWidth', 2);
    xlabel('S_0');
    ylabel('Price');
    title('Price of a Down-and-In Barrier Put Option vs S_0');
    legend('Down-and-In Barrier Put Option', 'Put Option');
    grid on;



end


