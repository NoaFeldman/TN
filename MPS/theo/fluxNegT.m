function res = fluxNegT(alpha, a, fixed)
    l = fixed(1);
    n = fixed(2);
    K = fixed(3);
    L = fixed(4);
    t = fixed(5);
    a = abs(a);
    epsilon = 1e-4; % a(1);
    [w11t, w22t, w33t, w12, w12t, w23, w23t, w13, w13t] = wDiffs(l, t, epsilon, L);
    dw1 = wDerivative(l, t, epsilon, L);    
    dw2 = wDerivative(0, t, epsilon, L);    
    dw3 = wDerivative(-l, t, epsilon, L);
    w11t = a(2) * w11t;
    w22t = a(2) * w22t;
    w33t = a(2) * w33t;
    w12  = a(2) * w12;
    w12t = a(2) * w12t;
    w23  = a(2) * w23;
    w23t = a(2) * w23t;
    w13  = a(2) * w13;
    w13t = a(2) * w13t;
    
    wpm = a(3);
    
    alpha1 = alpha / (2*pi);
    alpha2 = -2 .* alpha / (2*pi);
    alpha3 = alpha / (2*pi);
    order1 = (dw1).^(1 * K/n .* alpha1.^2) .* (dw2).^(1* K/n .* alpha2.^2) .* (dw3).^(1 * K/n .* alpha3.^2) .* ...
        (w11t).^(-1 * K/n .* alpha1.^2) .* (w22t).^(-1 * K/n .* alpha2.^2) .* (w33t).^(-1 * K/n .* alpha3.^2) .* ...
        (w12).^(2 * K/n .* alpha1 .* alpha2) .* (w12t).^(-2 * K/n .* alpha1 .* alpha2) .* ...
        (w23).^(2 * K/n .* alpha2 .* alpha3) .* (w23t).^(-2 * K/n .* alpha2 .* alpha3) .* ...
        (w13).^(2 * K/n .* alpha1 .* alpha3) .* (w13t).^(-2 * K/n .* alpha1 .* alpha3);
    order2 = 0;
    order2Shift = [0, 1, -1];
    for i = 1:3
        for j = 1:3
            for k = 1:3
                if i ~= j && i ~=k && k ~= j
                    alpha1Curr = alpha1 + order2Shift(i);
                    alpha2Curr = alpha2 + order2Shift(j);
                    alpha3Curr = alpha3 + order2Shift(k);
                    order2 = order2 + ...
                        (dw1).^(1 * K/n .* alpha1Curr.^2) .* (dw2).^(1 * K/n .* alpha2Curr.^2) .* (dw3).^(1 * K/n .* alpha3Curr.^2) .* ...
                        (w11t).^(-1 * K/n .* alpha1Curr.^2) .* (w22t).^(-1 * K/n .* alpha2Curr.^2) .* (w33t).^(-1 * K/n .* alpha3Curr.^2) .* ...
                        (w12).^(2 * K/n .* alpha1Curr .* alpha2Curr) .* (w12t).^(-2 * K/n .* alpha1Curr .* alpha2Curr) .* ...
                        (w23).^(2 * K/n .* alpha2Curr .* alpha3Curr) .* (w23t).^(-2 * K/n .* alpha2Curr .* alpha3Curr) .* ...
                        (w13).^(2 * K/n .* alpha1Curr .* alpha3Curr) .* (w13t).^(-2 * K/n .* alpha1Curr .* alpha3Curr);
                end
            end
        end
    end    
%     res = order1;
    res = order1 + wpm .* order2;
    zeroIndex = length(alpha) / 2 + 1/2;
%     res = (order1 + wpm .* order2) / (wpm * order2(zeroIndex) + order1(zeroIndex));
end