function res = fluxNegsT25(alphas, fixed)
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
    res = 0;
    for i = [1 3 4 6]
        alpha1Curr = alpha1;
        alpha2Curr = alpha2 + (i == 2) - (i == 5);
        alpha3Curr = alpha3;
        curr = ...
            (dw1).^(1 * K/n .* alpha1Curr.^2) .* (dw2).^(1 * K/n .* alpha2Curr.^2) .* (dw3).^(1 * K/n .* alpha3Curr.^2) .* ...
            (w11t).^(-1 * K/n .* alpha1Curr.^2) .* (w22t).^(-1 * K/n .* alpha2Curr.^2) .* (w33t).^(-1 * K/n .* alpha3Curr.^2) .* ...
            (w12).^(2 * K/n .* alpha1Curr .* alpha2Curr) .* (w12t).^(-2 * K/n .* alpha1Curr .* alpha2Curr) .* ...
            (w23).^(2 * K/n .* alpha2Curr .* alpha3Curr) .* (w23t).^(-2 * K/n .* alpha2Curr .* alpha3Curr) .* ...
            (w13).^(2 * K/n .* alpha1Curr .* alpha3Curr) .* (w13t).^(-2 * K/n .* alpha1Curr .* alpha3Curr); 
        % For alpha2, we consider V2 = V(-alpha)^2, and only
        % one of the operators gets the shift, so we get a
        % combinatoric factor of 2.
        res = res + 2 .* curr;
    end
end