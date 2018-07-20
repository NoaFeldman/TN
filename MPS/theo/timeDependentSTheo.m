function s = timeDependentSTheo(t, params, fixed)
    % stheo = -d(xS^(c(n - n^1)/6) N(0, sigma2) )/dn
    % xS is calculated from the sFull fit, and sigma2 from the P(NA) fit.
    L = fixed(1);
    NA = fixed(2);
    epsP = params(1);
    constP = params(2);
    epsS = params(3);
    c1 = params(4);
    sigma2 = getSigmaN(t, [epsP, constP], L);
    global c;
    c = 1;
%     xS = exp(-c1 * 12 / c) .* onePointFunc(t, epsS, L);
    xS = exp(-c1 * 6 / c) .* getOnePointFunc(t, epsS, L);
    s = stheo(xS, sigma2, NA);
end