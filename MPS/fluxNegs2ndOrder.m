function res = fluxNegs2ndOrder(alpha, a, fixed)
    res =  fluxNegT(alpha, [a(1) fixed(1) fixed(2) a(2)], fixed(3:6));
end