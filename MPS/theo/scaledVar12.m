function res = scaledVar12(t, a, fixed)
    l = fixed(1);
    n = fixed(2);
    epsilon = a(1);
    res = log (analyticallyContinuedExpr(l, t, epsilon).^(2 / n) .* ...
        sqrt(1 ./ eta12(l, t, epsilon).^(-4 /n))) + a(2);    
end