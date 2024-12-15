function Graph_Pricing(S0_values, V_Case1_Put, V_Case1_Perf, V_Case2_Put, V_Case2_Perf, proportionnal)
    figure;
    plot(S0_values, V_Case1_Put, 'b', 'LineWidth', 2);
    hold on;
    plot(S0_values, V_Case1_Perf, 'b--', 'LineWidth', 2);
    plot(S0_values, V_Case2_Put, 'r', 'LineWidth', 2);
    plot(S0_values, V_Case2_Perf, 'r--', 'LineWidth', 2);
    xlabel('S_0');
    ylabel('V(S_0, t_0)');
    if proportionnal == 1
        title('Option Prices vs Initial Stock Price S_0 with Proportionnal Barriers');
    else
        title('Option Prices vs Initial Stock Price S_0 with Fixed Barriers');
    end
    legend('Case 1 - Put Payoff', 'Case 1 - Performance Payoff', 'Case 2 - Put Payoff', 'Case 2 - Performance Payoff', 'Location', 'SouthEast');
    grid on;
end
