function res = fluxNeg2ndOrder(alpha, a, fixed)
    res = fluxNegT(alpha, [a 0], fixed);
end