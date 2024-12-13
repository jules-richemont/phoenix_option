% Function to calculate quantile using ordering algorithm
function x_alpha = Quantile(X, alpha)
    N_mc = length(X);
    X_sorted = sort(X);
    k = floor(N_mc * alpha);
    if k == 0
        x_alpha = X_sorted(1);
    else
        x_alpha = X_sorted(k);
    end
end