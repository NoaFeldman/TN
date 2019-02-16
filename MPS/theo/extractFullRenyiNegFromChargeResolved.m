function [rFullPred, rQ0Pred, rQ1Pred, rQ2Pred] = extractFullRenyiNegFromChargeResolved(t, l, epsilon, K, n, L, rQ2, rQ3)
    [w11t, w22t, w33t, w12, w12t, w23, w23t, w13, w13t] = wDiffs(l, t, epsilon, L);
    dw1 = wDerivative(l, t, epsilon, L);    
    dw2 = wDerivative(0, t, epsilon, L);    
    dw3 = wDerivative(-l, t, epsilon, L);
    leff = (dw1 ./ w11t .* dw3 ./ w33t .* dw2.^4 ./ w22t.^4 .* ...
        w12.^(-4) .* w12t.^(4) .* w23.^(-4) .* w23t.^(4) .* w13.^(2) .* w13t.^(-2)).^(-2);

    nuc = exp((-2*K .* log(rQ2./rQ3) / ((2^2 - 3^2) * n * pi^2)).^-1)./leff;
    rFullPred = rQ2 ./ (sqrt(pi * n ./ (2 * K .* log(leff.*nuc))) .* exp(-pi^2 * n * 2^2 ./ (2 * K .* log(leff.*nuc))));
    rQ0Pred = rFullPred .* sqrt(pi * n ./ (2 * K .* log(leff.*nuc)));
    rQ1Pred = rFullPred .* sqrt(pi * n ./ (2 * K .* log(leff.*nuc))).* exp(-pi^2 * n * 1^2 ./ (2 * K .* log(leff.*nuc)));
    rQ2Pred = rFullPred .* sqrt(pi * n ./ (2 * K .* log(leff.*nuc))).* exp(-pi^2 * n * 2^2 ./ (2 * K .* log(leff.*nuc)));
    
    rQ3Pred = rFullPred .* sqrt(pi * n ./ (2 * K .* log(leff.*nuc))).* exp(-pi^2 * n * 3^2 ./ (2 * K .* log(leff.*nuc)));
    rQ4Pred = rFullPred .* sqrt(pi * n ./ (2 * K .* log(leff.*nuc))).* exp(-pi^2 * n * 4^2 ./ (2 * K .* log(leff.*nuc)));
    rQ5Pred = rFullPred .* sqrt(pi * n ./ (2 * K .* log(leff.*nuc))).* exp(-pi^2 * n * 5^2 ./ (2 * K .* log(leff.*nuc)));
    
    rFullPred = rQ0Pred + 2 .* (rQ1Pred + rQ2 + rQ3 +rQ4Pred +rQ5Pred);
%     
%     nuc = exp((-2*K .* log(rQ3./rQ4) / ((3^2 - 4^2) * n * pi^2)).^-1)./leff;
%     rFullPred = rQ3 ./ (sqrt(pi * n ./ (2 * K .* log(leff.*nuc))) .* exp(-pi^2 * n * 3^2 ./ (2 * K .* log(leff.*nuc))));
%     rQ0Pred = rFullPred .* sqrt(pi * n ./ (2 * K .* log(leff.*nuc)));
%     rQ1Pred = rFullPred .* sqrt(pi * n ./ (2 * K .* log(leff.*nuc))).* exp(-pi^2 * n * 1^2 ./ (2 * K .* log(leff.*nuc)));
%     rQ2Pred = rFullPred .* sqrt(pi * n ./ (2 * K .* log(leff.*nuc))).* exp(-pi^2 * n * 2^2 ./ (2 * K .* log(leff.*nuc)));
%     
%     rQ3Pred = rFullPred .* sqrt(pi * n ./ (2 * K .* log(leff.*nuc))).* exp(-pi^2 * n * 3^2 ./ (2 * K .* log(leff.*nuc)));
%     rQ4Pred = rFullPred .* sqrt(pi * n ./ (2 * K .* log(leff.*nuc))).* exp(-pi^2 * n * 4^2 ./ (2 * K .* log(leff.*nuc)));
%     rQ5Pred = rFullPred .* sqrt(pi * n ./ (2 * K .* log(leff.*nuc))).* exp(-pi^2 * n * 5^2 ./ (2 * K .* log(leff.*nuc)));
%     
%     rFullPred = rQ0Pred + 2 .* (rQ1Pred + rQ2Pred + rQ3Pred +rQ4Pred +rQ5Pred);
end
