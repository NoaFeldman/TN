function res = fluxNegTUnfolded(alphas, a, fixed)
    l = fixed(1);
    n = fixed(2);
    K = fixed(3);
    L = fixed(4);
    t = fixed(5);
    a = abs(a);
    epsilon = 1; %  a(1);
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
    
    wpm1 = a(3);
    wpm2 = a(4);
    
    alpha1 = alphas / (2*pi);
    alpha2 = -2 .* alphas / (2*pi);
    alpha3 = alphas / (2*pi);
    order1 = (dw1).^(K/n .* alpha1.^2) .* (dw2).^(K/n .* alpha2.^2) .* (dw3).^(K/n .* alpha3.^2) .* ...
        (w11t).^(-1 * K/n .* alpha1.^2) .* (w22t).^(-1 * K/n .* alpha2.^2) .* (w33t).^(-1 * K/n .* alpha3.^2) .* ...
        (w12).^(2 * K/n .* alpha1 .* alpha2) .* (w12t).^(-2 * K/n .* alpha1 .* alpha2) .* ...
        (w23).^(2 * K/n .* alpha2 .* alpha3) .* (w23t).^(-2 * K/n .* alpha2 .* alpha3) .* ...
        (w13).^(2 * K/n .* alpha1 .* alpha3) .* (w13t).^(-2 * K/n .* alpha1 .* alpha3);
    order2 = 0;
    for i = 1:8
        for j = 1:8
            if j ~= i
                alpha1Curr = alpha1 + (i == 1) - (j == 1);
                alpha1CurrTilde = -alpha1 + (i == 2) - (j == 2);
                alpha2Curr = alpha2 + (i == 3) - (j == 3) + (i == 7) - (j == 7);
                alpha2CurrTilde = -alpha2 + (i == 4) - (j == 4) + (i == 8) - (j == 8);
                alpha3Curr = alpha3 + (i == 5) - (j == 5);
                alpha3CurrTilde = -alpha3 + (i == 6) - (j == 6);
                curr = ...
                    (dw1).^(K/n .* (alpha1Curr.^2 + alpha1CurrTilde.^2)) .* (dw2).^(K/n .* (alpha2Curr.^2 + alpha2CurrTilde.^2)) .* (dw3).^(K/n .* (alpha3Curr.^2 + alpha3CurrTilde.^2)) .* ...
                    (w11t).^(K/n .* alpha1Curr .* alpha1CurrTilde) .* ...
                    (w22t).^(K/n .* alpha2Curr .* alpha2CurrTilde) .* ...
                    (w33t).^(K/n .* alpha3Curr .* alpha3CurrTilde) .* ...
                    (w12).^(K/n .* alpha1Curr .* alpha2Curr) .* (w12).^(K/n .* alpha1CurrTilde .* alpha2CurrTilde) .* ...
                    (w12t).^(K/n .* alpha1Curr .* alpha2CurrTilde) .* (w12t).^(K/n .* alpha1CurrTilde .* alpha2Curr) .* ...
                    (w13).^(K/n .* alpha1Curr .* alpha3Curr) .* (w13).^(K/n .* alpha1CurrTilde .* alpha3CurrTilde) .* ...
                    (w13t).^(K/n .* alpha1Curr .* alpha3CurrTilde) .* (w13t).^(K/n .* alpha1CurrTilde .* alpha3Curr) .* ...
                    (w23).^(K/n .* alpha2Curr .* alpha3Curr) .* (w23).^(K/n .* alpha2CurrTilde .* alpha3CurrTilde) .* ...
                    (w23t).^(K/n .* alpha2Curr .* alpha3CurrTilde) .* (w23t).^(K/n .* alpha2CurrTilde .* alpha3Curr);                
                % For alpha2, we consider V2 = V(-alpha)^2, and only
                % one of the operators gets the shift, so we get a
                % combinatoric factor of 2.
                coeff = wpm1^2 * ((i == 3 || i == 4) + 1) * ((j == 3 || j == 4) + 1);
                order2 = order2 + coeff .* curr;
            end
        end
    end
    % Case 1 - shift 1 in four different places
    order3Case1 = 0;
    order3Case1Shift = [1, 1, -1, -1, 0, 0];
    inds = 1:6;
    for i = inds
        indsJ = inds(find(inds ~= i));
        for j = indsJ
            indsK = indsJ(find(indsJ ~= j));
            for k = indsK
                indsL = indsK(find(indsK ~= k));
                for l = indsL
                    indsM = indsL(find(indsL ~= l));
                    for m = indsM
                        indsO = indsM(find(indsM ~= m));
                        for o = indsO
                            alpha1Curr = alpha1 + order3Case1Shift(i);
                            alpha1CurrTilde = -alpha1 + order3Case1Shift(j);
                            alpha2Curr = alpha2 + order3Case1Shift(k);
                            alpha2CurrTilde = -alpha2 + order3Case1Shift(l);
                            alpha3Curr = alpha3 + order3Case1Shift(m);
                            alpha3CurrTilde = -alpha3 + order3Case1Shift(o);
                            curr = ...
                    (dw1).^(K/n .* (alpha1Curr.^2 + alpha1CurrTilde.^2)) .* (dw2).^(K/n .* (alpha2Curr.^2 + alpha2CurrTilde.^2)) .* (dw3).^(K/n .* (alpha3Curr.^2 + alpha3CurrTilde.^2)) .* ...
                                (w11t).^(K/n .* alpha1Curr .* alpha1CurrTilde) .* ...
                                (w22t).^(K/n .* alpha2Curr .* alpha2CurrTilde) .* ...
                                (w33t).^(K/n .* alpha3Curr .* alpha3CurrTilde) .* ...
                                (w12).^(K/n .* alpha1Curr .* alpha2Curr) .* (w12).^(K/n .* alpha1CurrTilde .* alpha2CurrTilde) .* ...
                                (w12t).^(K/n .* alpha1Curr .* alpha2CurrTilde) .* (w12t).^(K/n .* alpha1CurrTilde .* alpha2Curr) .* ...
                                (w13).^(K/n .* alpha1Curr .* alpha3Curr) .* (w13).^(K/n .* alpha1CurrTilde .* alpha3CurrTilde) .* ...
                                (w13t).^(K/n .* alpha1Curr .* alpha3CurrTilde) .* (w13t).^(K/n .* alpha1CurrTilde .* alpha3Curr) .* ...
                                (w23).^(K/n .* alpha2Curr .* alpha3Curr) .* (w23).^(K/n .* alpha2CurrTilde .* alpha3CurrTilde) .* ...
                                (w23t).^(K/n .* alpha2Curr .* alpha3CurrTilde) .* (w23t).^(K/n .* alpha2CurrTilde .* alpha3Curr);                
                            % For alpha2, we consider V2 = V(-alpha)^2, and only
                            % one of the operators gets the shift, so we get a
                            % combinatoric factor of 2.
                            combinatoricFactor = 1 * ((k < 5) + 1) * ((l < 5) + 1);
                            order3Case1 = order3Case1 + wpm1^4 .* combinatoricFactor .* curr;
                        end
                    end
                end
            end
        end
    end
    % Case 2 - two 1 shifts, one 2 shift.
    order3Case2 = 0;
    order3Case2Shift = [1, 1, -2, 0, 0, 0];
    inds = 1:6;
    for i = inds
        indsJ = inds(find(inds ~= i));
        for j = indsJ
            indsK = indsJ(find(indsJ ~= j));
            for k = indsK
                indsL = indsK(find(indsK ~= k));
                for l = indsL
                    indsM = indsL(find(indsL ~= l));
                    for m = indsM
                        indsO = indsM(find(indsM ~= m));
                        for o = indsO
                            alpha1Curr = alpha1 + order3Case2Shift(i);
                            alpha1CurrTilde = -alpha1 + order3Case2Shift(j);
                            alpha2Curr = alpha2 + order3Case2Shift(k);
                            alpha2CurrTilde = -alpha2 + order3Case2Shift(l);
                            alpha3Curr = alpha3 + order3Case2Shift(m);
                            alpha3CurrTilde = -alpha3 + order3Case2Shift(o);
                            curr = ...
                    (dw1).^(K/n .* (alpha1Curr.^2 + alpha1CurrTilde.^2)) .* (dw2).^(K/n .* (alpha2Curr.^2 + alpha2CurrTilde.^2)) .* (dw3).^(K/n .* (alpha3Curr.^2 + alpha3CurrTilde.^2)) .* ...
                                (w11t).^(K/n .* alpha1Curr .* alpha1CurrTilde) .* ...
                                (w22t).^(K/n .* alpha2Curr .* alpha2CurrTilde) .* ...
                                (w33t).^(K/n .* alpha3Curr .* alpha3CurrTilde) .* ...
                                (w12).^(K/n .* alpha1Curr .* alpha2Curr) .* (w12).^(K/n .* alpha1CurrTilde .* alpha2CurrTilde) .* ...
                                (w12t).^(K/n .* alpha1Curr .* alpha2CurrTilde) .* (w12t).^(K/n .* alpha1CurrTilde .* alpha2Curr) .* ...
                                (w13).^(K/n .* alpha1Curr .* alpha3Curr) .* (w13).^(K/n .* alpha1CurrTilde .* alpha3CurrTilde) .* ...
                                (w13t).^(K/n .* alpha1Curr .* alpha3CurrTilde) .* (w13t).^(K/n .* alpha1CurrTilde .* alpha3Curr) .* ...
                                (w23).^(K/n .* alpha2Curr .* alpha3Curr) .* (w23).^(K/n .* alpha2CurrTilde .* alpha3CurrTilde) .* ...
                                (w23t).^(K/n .* alpha2Curr .* alpha3CurrTilde) .* (w23t).^(K/n .* alpha2CurrTilde .* alpha3Curr);                
                            % For alpha2, we consider V2 = V(-alpha)^2, and only
                            % one of the operators gets the shift, so we get a
                            % combinatoric factor of 2.
                            coeff = 1;
                            if (k == 3 || l == 3)
                                coeff = 2 * wpm1^2 * wpm2 + wpm1^4;
                                if (k < 3 || l < 3)
                                    coeff = coeff*2;
                                end
                            elseif (k > 3 && l > 3)
                                coeff = wpm1^2 * wpm2;
                            elseif (k < 3 || l < 3)
                                coeff = 2 * wpm1^2 * wpm2;
                            else
                                pause(1);
                            end
                            order3Case2 = order3Case2 +  coeff .* curr;
                            
                            alpha1Curr = alpha1 - order3Case2Shift(i);
                            alpha1CurrTilde = -alpha1 - order3Case2Shift(j);
                            alpha2Curr = alpha2 - order3Case2Shift(k);
                            alpha2CurrTilde = -alpha2 - order3Case2Shift(l);
                            alpha3Curr = alpha3 - order3Case2Shift(m);
                            alpha3CurrTilde = -alpha3 - order3Case2Shift(o);
                            curr = ...
                    (dw1).^(K/n .* (alpha1Curr.^2 + alpha1CurrTilde.^2)) .* (dw2).^(K/n .* (alpha2Curr.^2 + alpha2CurrTilde.^2)) .* (dw3).^(K/n .* (alpha3Curr.^2 + alpha3CurrTilde.^2)) .* ...
                                (w11t).^(K/n .* alpha1Curr .* alpha1CurrTilde) .* ...
                                (w22t).^(K/n .* alpha2Curr .* alpha2CurrTilde) .* ...
                                (w33t).^(K/n .* alpha3Curr .* alpha3CurrTilde) .* ...
                                (w12).^(K/n .* alpha1Curr .* alpha2Curr) .* (w12).^(K/n .* alpha1CurrTilde .* alpha2CurrTilde) .* ...
                                (w12t).^(K/n .* alpha1Curr .* alpha2CurrTilde) .* (w12t).^(K/n .* alpha1CurrTilde .* alpha2Curr) .* ...
                                (w13).^(K/n .* alpha1Curr .* alpha3Curr) .* (w13).^(K/n .* alpha1CurrTilde .* alpha3CurrTilde) .* ...
                                (w13t).^(K/n .* alpha1Curr .* alpha3CurrTilde) .* (w13t).^(K/n .* alpha1CurrTilde .* alpha3Curr) .* ...
                                (w23).^(K/n .* alpha2Curr .* alpha3Curr) .* (w23).^(K/n .* alpha2CurrTilde .* alpha3CurrTilde) .* ...
                                (w23t).^(K/n .* alpha2Curr .* alpha3CurrTilde) .* (w23t).^(K/n .* alpha2CurrTilde .* alpha3Curr);                
                            % For alpha2, we consider V2 = V(-alpha)^2, and only
                            % one of the operators gets the shift, so we get a
                            % combinatoric factor of 2.
                            coeff = 1;
                            if (k == 3 || l == 3)
                                coeff = 2 * wpm1^2 * wpm2 + wpm1^4;
                                if (k < 3 || l < 3)
                                    coeff = coeff*2;
                                end
                            elseif (k > 3 && l > 3)
                                coeff = wpm1^2 * wpm2;
                            elseif (k < 3 || l < 3)
                                coeff = 2 * wpm1^2 * wpm2;
                            else
                                pause(1);
                            end
                            order3Case2 = order3Case2 +  coeff .* curr;
                        end
                    end
                end
            end
        end
    end
    % Case 3 - two 2 shifts.
    order3Case3 = 0;
    order3Case3Shift = [2, -2, 0, 0, 0, 0];
    inds = 1:6;
    for i = inds
        indsJ = inds(find(inds ~= i));
        for j = indsJ
            indsK = indsJ(find(indsJ ~= j));
            for k = indsK
                indsL = indsK(find(indsK ~= k));
                for l = indsL
                    indsM = indsL(find(indsL ~= l));
                    for m = indsM
                        indsO = indsM(find(indsM ~= m));
                        for o = indsO
                            alpha1Curr = alpha1 + order3Case3Shift(i);
                            alpha1CurrTilde = -alpha1 + order3Case3Shift(j);
                            alpha2Curr = alpha2 + order3Case3Shift(k);
                            alpha2CurrTilde = -alpha2 + order3Case3Shift(l);
                            alpha3Curr = alpha3 + order3Case3Shift(m);
                            alpha3CurrTilde = -alpha3 + order3Case3Shift(o);
                            curr = ...
                    (dw1).^(K/n .* (alpha1Curr.^2 + alpha1CurrTilde.^2)) .* (dw2).^(K/n .* (alpha2Curr.^2 + alpha2CurrTilde.^2)) .* (dw3).^(K/n .* (alpha3Curr.^2 + alpha3CurrTilde.^2)) .* ...
                                (w11t).^(K/n .* alpha1Curr .* alpha1CurrTilde) .* ...
                                (w22t).^(K/n .* alpha2Curr .* alpha2CurrTilde) .* ...
                                (w33t).^(K/n .* alpha3Curr .* alpha3CurrTilde) .* ...
                                (w12).^(K/n .* alpha1Curr .* alpha2Curr) .* (w12).^(K/n .* alpha1CurrTilde .* alpha2CurrTilde) .* ...
                                (w12t).^(K/n .* alpha1Curr .* alpha2CurrTilde) .* (w12t).^(K/n .* alpha1CurrTilde .* alpha2Curr) .* ...
                                (w13).^(K/n .* alpha1Curr .* alpha3Curr) .* (w13).^(K/n .* alpha1CurrTilde .* alpha3CurrTilde) .* ...
                                (w13t).^(K/n .* alpha1Curr .* alpha3CurrTilde) .* (w13t).^(K/n .* alpha1CurrTilde .* alpha3Curr) .* ...
                                (w23).^(K/n .* alpha2Curr .* alpha3Curr) .* (w23).^(K/n .* alpha2CurrTilde .* alpha3CurrTilde) .* ...
                                (w23t).^(K/n .* alpha2Curr .* alpha3CurrTilde) .* (w23t).^(K/n .* alpha2CurrTilde .* alpha3Curr);                
                            % For alpha2, we consider V2 = V(-alpha)^2, and only
                            % one of the operators gets the shift, so we get a
                            % combinatoric factor of 2.
                            if (k > 2 && l > 2)
                                coeff = wpm2^2;
                            elseif (k > 2 || l > 2)
                                coeff = 2 * wpm2^2 + wpm2 * wpm1^2;
                            else
                                coeff = 4 * wpm2^2 + 4 * wpm2 * wpm1^2 + wpm1^4;
                            end
                            order3Case3 = order3Case3 +  coeff .* curr;
                        end
                    end
                end
            end
        end
    end
    
    zeroIndex = length(alphas) / 2 + 1/2;
    res = (order1 + order2 + order3Case1 + order3Case2 + order3Case3) / ...
        (order1(zeroIndex) + order2(zeroIndex) + order3Case1(zeroIndex) + order3Case2(zeroIndex) + order3Case3(zeroIndex));
end