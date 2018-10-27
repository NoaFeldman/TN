function res = renyiNegT(alpha, a, fixed)
    l = fixed(1);
    t = fixed(2);
    n = fixed(3);
    K = fixed(4);
    epsilon = a(1);
    [w11t, w22t, w33t, w12, w12t, w23, w23t, w13, w13t] = getWDiffs(l, t, epsilon);
    dw1 = wDerivative(l, t, epsilon);    
    dw2 = wDerivative(0, t, epsilon);    
    dw3 = wDerivative(-l, t, epsilon);
    alpha1 = alpha / (2*pi);
    alpha2 = -2 .* alpha / (2*pi);
    alpha3 = alpha / (2*pi);
    order1 = (dw1).^(K/n .* alpha1.^2) .* (dw2).^(K/n .* alpha2.^2) .* (dw3).^(K/n .* alpha3.^2) .* ...
        (w11t).^(-K/n .* alpha1.^2) .* (w22t).^(-K/n .* alpha2.^2) .* (w33t).^(-K/n .* alpha3.^2) .* ...
        (w12).^(K/n .* alpha1 .* alpha2) .* (w12t).^(-K/n .* alpha1 .* alpha2) .* ...
        (w23).^(K/n .* alpha2 .* alpha3) .* (w23t).^(-K/n .* alpha2 .* alpha3) .* ...
        (w13).^(K/n .* alpha1 .* alpha3) .* (w13t).^(-K/n .* alpha1 .* alpha3);
    res = order1;
end