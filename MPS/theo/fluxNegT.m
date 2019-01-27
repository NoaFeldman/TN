function res = fluxNegT(alphas, a, fixed)
% i = 30, plot(alphas, fluxNegT(alphas, [1 1.55 0.6635 0.9 0.2], [l n K L t(i)]));
    l = fixed(1);
    n = fixed(2);
    K = fixed(3);
    L = fixed(4);
    t = fixed(5);
    a = abs(a);
    epsilon = 1; % a(1);
    nuc = a(2);
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
    
    wpm2 = a(3);
    wpm3 = a(4);
    
    alpha1 = alphas / (2*pi);
    alpha2 = -2 .* alphas / (2*pi);
    alpha3 = alphas / (2*pi);
    order1 = (dw1).^(1 * K/n .* alpha1.^2) .* (dw2).^(1* K/n .* alpha2.^2) .* (dw3).^(1 * K/n .* alpha3.^2) .* ...
        (w11t).^(-1 * K/n .* alpha1.^2) .* (w22t).^(-1 * K/n .* alpha2.^2) .* (w33t).^(-1 * K/n .* alpha3.^2) .* ...
        (w12).^(2 * K/n .* alpha1 .* alpha2) .* (w12t).^(-2 * K/n .* alpha1 .* alpha2) .* ...
        (w23).^(2 * K/n .* alpha2 .* alpha3) .* (w23t).^(-2 * K/n .* alpha2 .* alpha3) .* ...
        (w13).^(2 * K/n .* alpha1 .* alpha3) .* (w13t).^(-2 * K/n .* alpha1 .* alpha3);
    order2 = 0;
    for i = 1:6
        alpha1Curr = alpha1 + (i == 1) - (i == 4);
        alpha2Curr = alpha2 + (i == 2) - (i == 5);
        alpha3Curr = alpha3 + (i == 3) - (i == 6);
        curr = ...
            (dw1).^(1 * K/n .* alpha1Curr.^2) .* (dw2).^(1 * K/n .* alpha2Curr.^2) .* (dw3).^(1 * K/n .* alpha3Curr.^2) .* ...
            (w11t).^(-1 * K/n .* alpha1Curr.^2) .* (w22t).^(-1 * K/n .* alpha2Curr.^2) .* (w33t).^(-1 * K/n .* alpha3Curr.^2) .* ...
            (w12).^(2 * K/n .* alpha1Curr .* alpha2Curr) .* (w12t).^(-2 * K/n .* alpha1Curr .* alpha2Curr) .* ...
            (w23).^(2 * K/n .* alpha2Curr .* alpha3Curr) .* (w23t).^(-2 * K/n .* alpha2Curr .* alpha3Curr) .* ...
            (w13).^(2 * K/n .* alpha1Curr .* alpha3Curr) .* (w13t).^(-2 * K/n .* alpha1Curr .* alpha3Curr); 
        % For alpha2, we consider V2 = V(-alpha)^2, and only
        % one of the operators gets the shift, so we get a
        % combinatoric factor of 2.
        if i == 2 || i == 5
            order2 = order2 + 2 * wpm2 .* curr;
        else
            order2 = order2 + wpm2 .* curr;
        end
    end
    
    order3 = 0;
    for i = 1:6
        for j = 1:6
            if (i ~= j && abs(i - j) ~= 3)
                alpha1Curr = alpha1 + (i == 1) - (i == 4) + (j == 1) - (j == 4);
                alpha2Curr = alpha2 + (i == 2) - (i == 5) + (j == 2) - (j == 5);
                alpha3Curr = alpha3 + (i == 3) - (i == 6) + (j == 3) - (j == 6);
                curr = ...
                    (dw1).^(1 * K/n .* alpha1Curr.^2) .* (dw2).^(1 * K/n .* alpha2Curr.^2) .* (dw3).^(1 * K/n .* alpha3Curr.^2) .* ...
                    (w11t).^(-1 * K/n .* alpha1Curr.^2) .* (w22t).^(-1 * K/n .* alpha2Curr.^2) .* (w33t).^(-1 * K/n .* alpha3Curr.^2) .* ...
                    (w12).^(2 * K/n .* alpha1Curr .* alpha2Curr) .* (w12t).^(-2 * K/n .* alpha1Curr .* alpha2Curr) .* ...
                    (w23).^(2 * K/n .* alpha2Curr .* alpha3Curr) .* (w23t).^(-2 * K/n .* alpha2Curr .* alpha3Curr) .* ...
                    (w13).^(2 * K/n .* alpha1Curr .* alpha3Curr) .* (w13t).^(-2 * K/n .* alpha1Curr .* alpha3Curr); 
                % For alpha2, we consider V2 = V(-alpha)^2, and only
                % one of the operators gets the shift, so we get a
                % combinatoric factor of 2.
                combFactor = 1;
                if i == 2 || i == 5
                    combFactor = combFactor * 2;
                end
                if j == 2 || j == 5
                    combFactor = combFactor * 2;                
                end
                order3 = order3 + combFactor .* wpm2^2 .* curr;
            end
        end
    end
    
    for i = 1:6
        alpha1Curr = alpha1 + 2*((i == 1) - (i == 4));
        alpha2Curr = alpha2 + 2*((i == 2) - (i == 5));
        alpha3Curr = alpha3 + 2*((i == 3) - (i == 6));
        curr = ...
            (dw1).^(1 * K/n .* alpha1Curr.^2) .* (dw2).^(1 * K/n .* alpha2Curr.^2) .* (dw3).^(1 * K/n .* alpha3Curr.^2) .* ...
            (w11t).^(-1 * K/n .* alpha1Curr.^2) .* (w22t).^(-1 * K/n .* alpha2Curr.^2) .* (w33t).^(-1 * K/n .* alpha3Curr.^2) .* ...
            (w12).^(2 * K/n .* alpha1Curr .* alpha2Curr) .* (w12t).^(-2 * K/n .* alpha1Curr .* alpha2Curr) .* ...
            (w23).^(2 * K/n .* alpha2Curr .* alpha3Curr) .* (w23t).^(-2 * K/n .* alpha2Curr .* alpha3Curr) .* ...
            (w13).^(2 * K/n .* alpha1Curr .* alpha3Curr) .* (w13t).^(-2 * K/n .* alpha1Curr .* alpha3Curr); 
        % For alpha2, we consider V2 = V(-alpha)^2, and only
        % one of the operators gets the shift, so we get a
        % combinatoric factor of 2.
        combFactor = 1;
        if i == 2 || i == 5
            combFactor = combFactor * 2;
        end
        order3 = order3 + combFactor .* wpm3 .* curr;
    end
    
    order4 = 0;
    for i = 1:6
        for j = 1:6
            for k = 1:6
                if (i ~= j && abs(i - j) ~= 3 && i ~= k && abs(i - k) ~= 3 && k ~= j && abs(k - j) ~= 3)
                    alpha1Curr = alpha1 + (i == 1) - (i == 4) + (j == 1) - (j == 4) + (k == 1) - (k == 4);
                    alpha2Curr = alpha2 + (i == 2) - (i == 5) + (j == 2) - (j == 5) + (k == 2) - (k == 5);
                    alpha3Curr = alpha3 + (i == 3) - (i == 6) + (j == 3) - (j == 6) + (k == 3) - (k == 6);
                    curr = ...
                        (dw1).^(1 * K/n .* alpha1Curr.^2) .* (dw2).^(1 * K/n .* alpha2Curr.^2) .* (dw3).^(1 * K/n .* alpha3Curr.^2) .* ...
                        (w11t).^(-1 * K/n .* alpha1Curr.^2) .* (w22t).^(-1 * K/n .* alpha2Curr.^2) .* (w33t).^(-1 * K/n .* alpha3Curr.^2) .* ...
                        (w12).^(2 * K/n .* alpha1Curr .* alpha2Curr) .* (w12t).^(-2 * K/n .* alpha1Curr .* alpha2Curr) .* ...
                        (w23).^(2 * K/n .* alpha2Curr .* alpha3Curr) .* (w23t).^(-2 * K/n .* alpha2Curr .* alpha3Curr) .* ...
                        (w13).^(2 * K/n .* alpha1Curr .* alpha3Curr) .* (w13t).^(-2 * K/n .* alpha1Curr .* alpha3Curr); 
                    % For alpha2, we consider V2 = V(-alpha)^2, and only
                    % one of the operators gets the shift, so we get a
                    % combinatoric factor of 2.
                    combFactor = 1;
                    if i == 2 || i == 5
                        combFactor = combFactor * 2;
                    end
                    if j == 2 || j == 5
                        combFactor = combFactor * 2;                
                    end
                    if k == 2 || k == 5
                        combFactor = combFactor * 2;                
                    end
                    order4 = order4 + combFactor .* wpm2^3 .* curr;
                end
            end
        end
    end
    
    for i = 1:6
        for j = 1:6
            if (i ~= j && abs(i - j) ~= 3)
                alpha1Curr = alpha1 + (i == 1) - (i == 4) + 2*((j == 1) - (j == 4));
                alpha2Curr = alpha2 + (i == 2) - (i == 5) + 2*((j == 2) - (j == 5));
                alpha3Curr = alpha3 + (i == 3) - (i == 6) + 2*((j == 3) - (j == 6));
                curr = ...
                    (dw1).^(1 * K/n .* alpha1Curr.^2) .* (dw2).^(1 * K/n .* alpha2Curr.^2) .* (dw3).^(1 * K/n .* alpha3Curr.^2) .* ...
                    (w11t).^(-1 * K/n .* alpha1Curr.^2) .* (w22t).^(-1 * K/n .* alpha2Curr.^2) .* (w33t).^(-1 * K/n .* alpha3Curr.^2) .* ...
                    (w12).^(2 * K/n .* alpha1Curr .* alpha2Curr) .* (w12t).^(-2 * K/n .* alpha1Curr .* alpha2Curr) .* ...
                    (w23).^(2 * K/n .* alpha2Curr .* alpha3Curr) .* (w23t).^(-2 * K/n .* alpha2Curr .* alpha3Curr) .* ...
                    (w13).^(2 * K/n .* alpha1Curr .* alpha3Curr) .* (w13t).^(-2 * K/n .* alpha1Curr .* alpha3Curr); 
                % For alpha2, we consider V2 = V(-alpha)^2, and only
                % one of the operators gets the shift, so we get a
                % combinatoric factor of 2.
                combFactor = 1;
                if i == 2 || i == 5
                    combFactor = combFactor * 2;
                end
                if j == 2 || j == 5
                    combFactor = combFactor * 2;                
                end
                order4 = order4 + combFactor .* wpm2*wpm3 .* curr;
            end
        end
    end
    
%     res = order1 + wpm^2 .* order2;
    zeroIndex = length(alphas) / 2 + 1/2;
%     res = (order1(zeroIndex) + order2(zeroIndex) + order3(zeroIndex) + order4(zeroIndex));
    res = (order1 + order2 + order3 + order4) / ...
        (order1(zeroIndex) + order2(zeroIndex) + order3(zeroIndex) + order4(zeroIndex));
end