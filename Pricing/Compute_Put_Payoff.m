function Payoff = Compute_Put_Payoff(S, K, Pi0, C_Ph, C_Y, T, r, B_Ph, B_Y, B_Put, delta_t)
    timesteps = length(S) - 1;
    t = delta_t * (1:timesteps);
    PV = 0; % Présent Value
    autocall = false;

    % Boucle sur les barrières intermédiaires
    for i = 1:timesteps
        if S(i) > B_Ph
            % Autocall déclenché
            PV = (Pi0 + C_Ph) * exp(-r * t(i));
            autocall = true;
            break;
        elseif S(i) >= B_Y && S(i) <= B_Ph
            % Paiement des coupons
            PV = PV + C_Y * exp(-r * t(i));
        end
    end

    if ~autocall
        % Payoff à maturité
        S_N = S(end);
        if S_N > B_Ph
            PV = PV + (Pi0 + C_Ph) * exp(-r * T);
        elseif S_N >= B_Y && S_N <= B_Ph
            PV = PV + (Pi0 + C_Y) * exp(-r * T);
        elseif S_N > B_Put && S_N < B_Y
            PV = PV + Pi0 * exp(-r * T);
        elseif S_N <= B_Put
            % Payoff Put
            Payoff = max((K - S_N) / S(1), 0);
            PV = PV + Payoff * exp(-r * T);
        end
    end

    Payoff = PV; % Renvoyer la valeur finale actualisée
end
