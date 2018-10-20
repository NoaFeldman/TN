function res = renyiFluxSigma2(t, a, fixed)
    l = fixed(1);
    K = fixed(2);
    n = fixed(3);
    epsilon = a(1);
    res = log (analyticallyContinuedExpr(l, t, epsilon).^(2 * K / n) .* analyticallyContinuedExpr(-l, t, epsilon).^(2 * K / n) .* ...
        analyticallyContinuedExpr(0, t, epsilon).^(-4 * K / n) .* ...
        sqrt(eta13(l, t, epsilon).^(-8 * K / n) ./ eta12(l, t, epsilon).^(-8 * K /n))) + a(2);    
end