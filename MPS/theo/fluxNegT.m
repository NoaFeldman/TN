function res = fluxNegT(alpha, a, fixed)
    l = fixed(1);
    n = fixed(2);
    K = fixed(3);
    t = fixed(4);
    a = real(a);
    epsilon = 30; %a(1);
    [w11t, w22t, w33t, w12, w12t, w23, w23t, w13, w13t] = wDiffs(l/2, t, epsilon);
    dw1 = wDerivative(l/2, t, epsilon);    
    dw2 = wDerivative(0, t, epsilon);    
    dw3 = wDerivative(-l/2, t, epsilon);
    a(2) = abs(a(2));
    a(3) = abs(a(3)); % abs(a(2)); %
    w11t = a(2)^2 * w11t;
    w22t = a(3)^2 * w22t;
    w33t = a(2)^2 * w33t;
    w12 = a(2)*a(3) * w12;
    w12t = a(2)*a(3) * w12t;
    w23 = a(2)*a(3) * w23;
    w23t = a(2)*a(3) * w23t;
    w13 = a(2)^2 * w13;
    w13t = a(2)^2 * w13t;
    dw1 = a(2)^2 * dw1;
    dw2 = a(3)^2 * dw2;
    dw3 = a(2)^2 * dw3;
    wpm = a(4);
    
    alpha1 = alpha / (2*pi);
    alpha2 = -2 .* alpha / (2*pi);
    alpha3 = alpha / (2*pi);
    order1 = (dw1).^(4 * 1/2 * K/n .* alpha1.^2) .* (dw2).^(4 * 1/2 * K/n .* alpha2.^2) .* (dw3).^(4 * 1/2 * K/n .* alpha3.^2) .* ...
        (w11t).^(-4 * K/n .* alpha1.^2) .* (w22t).^(-4 * K/n .* alpha2.^2) .* (w33t).^(-4 * K/n .* alpha3.^2) .* ...
        (w12).^(4 * K/n .* alpha1 .* alpha2) .* (w12t).^(-4 * K/n .* alpha1 .* alpha2) .* ...
        (w23).^(4 * K/n .* alpha2 .* alpha3) .* (w23t).^(-4 * K/n .* alpha2 .* alpha3) .* ...
        (w13).^(4 * K/n .* alpha1 .* alpha3) .* (w13t).^(-4 * K/n .* alpha1 .* alpha3);
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
                        (dw1).^(4/2 * K/n .* alpha1Curr.^2) .* (dw2).^(4/2 * K/n .* alpha2Curr.^2) .* (dw3).^(4/2 * K/n .* alpha3Curr.^2) .* ...
        (w11t).^(-4 * K/n .* alpha1Curr.^2) .* (w22t).^(-4 * K/n .* alpha2Curr.^2) .* (w33t).^(-4 * K/n .* alpha3Curr.^2) .* ...
        (w12).^(4 * K/n .* alpha1Curr .* alpha2Curr) .* (w12t).^(-4 * K/n .* alpha1Curr .* alpha2Curr) .* ...
        (w23).^(4 * K/n .* alpha2Curr .* alpha3Curr) .* (w23t).^(-4 * K/n .* alpha2Curr .* alpha3Curr) .* ...
        (w13).^(4 * K/n .* alpha1Curr .* alpha3Curr) .* (w13t).^(-4 * K/n .* alpha1Curr .* alpha3Curr);
                end
            end
        end
    end    
%     res = order1;
%     res = order1 + wpm .* order2;
    zeroIndex = length(alpha) / 2 + 1/2;;
    res = (order1 + wpm .* order2) / (wpm * order2(zeroIndex) + order1(zeroIndex));
end