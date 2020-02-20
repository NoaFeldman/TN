function res = fluxNegsT0(alphas, fixed)
    l = fixed(1);
    n = fixed(2);
    K = fixed(3);
    L = fixed(4);
    t = fixed(5);
    epsilon = 1; % a(1);
    nuc = 1;
    [w11t, w22t, w33t, w12, w12t, w23, w23t, w13, w13t] = wDiffs(l, t, epsilon, L);
    dw1 = wDerivative(l, t, epsilon, L);    
    dw2 = wDerivative(0, t, epsilon, L);    
    dw3 = wDerivative(-l, t, epsilon, L);
    w11t = nuc * w11t;
    w22t = nuc * w22t;
    w33t = nuc * w33t;
    w12  = nuc * w12;
    w12t = nuc * w12t;
    w23  = nuc * w23;
    w23t = nuc * w23t;
    w13  = nuc * w13;
    w13t = nuc * w13t;
    
    alpha1 = alphas / (2*pi);
    alpha2 = -2 .* alphas / (2*pi);
    alpha3 = alphas / (2*pi);
    res = (dw1).^(1 * K/n .* alpha1.^2) .* (dw2).^(1* K/n .* alpha2.^2) .* (dw3).^(1 * K/n .* alpha3.^2) .* ...
        (w11t).^(-1 * K/n .* alpha1.^2) .* (w22t).^(-1 * K/n .* alpha2.^2) .* (w33t).^(-1 * K/n .* alpha3.^2) .* ...
        (w12).^(2 * K/n .* alpha1 .* alpha2) .* (w12t).^(-2 * K/n .* alpha1 .* alpha2) .* ...
        (w23).^(2 * K/n .* alpha2 .* alpha3) .* (w23t).^(-2 * K/n .* alpha2 .* alpha3) .* ...
        (w13).^(2 * K/n .* alpha1 .* alpha3) .* (w13t).^(-2 * K/n .* alpha1 .* alpha3);
end