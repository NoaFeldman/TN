function res = fluxNeg1stOrder(alpha, a, fixed)
    res = fluxNegT(alpha, [a 0], fixed);
end