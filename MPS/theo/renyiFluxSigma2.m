function res = renyiFluxSigma2(t, a, fixed)
    l = fixed(1);
    n = fixed(2);
    epsilon = a(1);
    res = log (analyticallyContinuedExpr(l, t, epsilon).^(2 / n) .* analyticallyContinuedExpr(-l, t, epsilon).^(2 / n) .* ...
        analyticallyContinuedExpr(0, t, epsilon).^(-4 / n) .* ...
        sqrt(eta13(l, t, epsilon).^(-8 / n) ./ eta12(l, t, epsilon).^(-8 /n))) + a(2);    
end