function Payoff = Compute_Payoff(S, K, Pi0, C_Ph, C_Y, T, r, B_Ph, B_Y, B_Put, delta_t, put_or_perf)
    % Initialisation
    timesteps = length(S) - 1; % Nombre de points moins le point initial
    t = delta_t * (0:timesteps) / timesteps * T; % Correspondance des pas de temps
    PV = 0; % Valeur actuelle des flux
    autocall = false;

    % Déterminer les indices des dates d'observation
    observation_date = round(linspace(1, length(S), T / delta_t + 1)); % Espacement uniforme

    % Boucle sur les dates d'observation
    for i = observation_date
        fprintf('i: %d\n', i);
        if S(i) > B_Ph
            % Autocall déclenché
            PV = (Pi0 + C_Ph) * exp(-r * t(i));
            autocall = true;
            fprintf('time: %d\n', t(i));
            break;
        elseif S(i) >= B_Y && S(i) <= B_Ph
            % Paiement d’un coupon
            PV = PV + C_Y * exp(-r * t(i));
            fprintf('time: %d\n', t(i));
        end
    end

    if ~autocall % If the autocall is not triggered, ie the autocall barrier is not crossed
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
            if put_or_perf == 1
                Payoff = max((K - S_N) / S(1), 0);
            else
                Payoff = S_N / S(1);
            end
            PV = PV + Payoff * exp(-r * T);
        end
    end

    Payoff = PV; % Renvoyer la valeur finale actualisée
end
