function Payoff = Compute_Perf_Payoff(S, Pi0, C_Ph, C_Y, T, r, B_Ph, B_Y, B_Put, delta_t)
    timesteps = length(S) - 1;
    t = delta_t * (1:timesteps);
    PV = 0;
    autocall = false;

    for i = 1:timesteps
        if S(i) > B_Ph
            PV = (Pi0 + C_Ph) * exp(-r * t(i));
            autocall = true;
            break;
        elseif S(i) >= B_Y && S(i) <= B_Ph
            PV = PV + C_Y * exp(-r * t(i));
        end
    end

    if ~autocall
        S_N = S(end);
        if S_N <= B_Put
            Payoff = (S_N / S(1)) * exp(-r * T);
        else
            Payoff = PV + Pi0 * exp(-r * T);
        end
    else
        Payoff = PV;
    end
end
