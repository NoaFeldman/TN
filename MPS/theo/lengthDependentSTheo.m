function s = lengthDependentSTheo(Ls, params, NA)
    X_scaled = exp(params(1)*6)./Ls;
    sigma2 = log(Ls)*params(2) + params(3);
    s = stheo(X_scaled, sigma2, NA);
end