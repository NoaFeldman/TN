function res = scaledVar13(t, a, fixed)
    l = fixed(1);
    n = fixed(2);
    epsilon = a(1);
    res = analyticallyContinuedExpr(0, t, epsilon).^(-4 / n) .* ...
        sqrt(eta13(l, t, epsilon).^(-8 / n));    
end