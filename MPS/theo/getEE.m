function S = getEE(f, x, L)
    dn = 1e-4;
    snTop = getSNA(1 + dn / 2, f, x, L);
    snBottom = getSNA(1 - dn / 2, f, x, L);
    S = -(snTop - snBottom) / dn;
end