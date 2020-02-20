function res = fluxNegTPeriodic(alphas, a, fixed, periodNum)
    res = fluxNegT(alphas, a, fixed);
    for i = 1:periodNum
        res = res + fluxNegT(alphas + i * 2*pi, a, fixed);
        res = res + fluxNegT(alphas - i * 2*pi, a, fixed);
    end
    zeroIndex = find(alphas == 0);
end
    