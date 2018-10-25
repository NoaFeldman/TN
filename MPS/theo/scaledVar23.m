function res = scaledVar23(t, a, fixed)
    l = fixed(1);
    n = fixed(2);
    epsilon = a(1);
    res = analyticallyContinuedExpr(-l, t, epsilon).^(2 / n) .* ...
        sqrt(eta13(l, t, epsilon).^(-2 / n) ./ eta12(l, t, epsilon).^(-4 /n));    
end