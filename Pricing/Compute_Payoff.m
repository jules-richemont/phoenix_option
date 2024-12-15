function Payoff = Compute_Payoff(S,K,Pi0,C_Ph,C_Y,T,r,B_Ph,B_Y,B_Put,delta_t,put_or_perf)
    % Initialization
    timesteps = length(S)-1; % Number of points minus the initial point
    t = delta_t*(0:timesteps)/timesteps*T; % Time step correspondence
    PV = 0; % Present value of cash flows
    autocall = false;

    % Determine the indices of the observation dates
    observation_date = round(linspace(1,length(S),T/delta_t+1)); % Uniform spacing

    % Loop over the observation dates
    for i = observation_date
        if S(i)>B_Ph
            % Autocall triggered
            PV = (Pi0+C_Ph)*exp(-r*t(i));
            autocall = true;
            break;
        elseif S(i)>=B_Y && S(i)<=B_Ph
            % Coupon payment
            PV = PV + C_Y*exp(-r*t(i));
        end
    end

    if ~autocall % If the autocall is not triggered, i.e., the autocall barrier is not crossed
        % Payoff at maturity
        S_N = S(end);
        if S_N>B_Ph
            PV = PV + (Pi0+C_Ph)*exp(-r*T);
        elseif S_N>=B_Y && S_N<=B_Ph
            PV = PV + (Pi0+C_Y)*exp(-r*T);
        elseif S_N>B_Put && S_N<B_Y
            PV = PV + Pi0*exp(-r*T);
        elseif S_N<=B_Put
            % Put Payoff
            if put_or_perf == 1
                Payoff = max((K-S_N)/S(1),0);
            else
                Payoff = S_N/S(1);
            end
            PV = PV + Payoff*exp(-r*T);
        end
    end

    Payoff = PV; % Return the final discounted value
end
