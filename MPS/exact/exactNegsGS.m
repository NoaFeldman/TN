function exactNegsGS(L, l, filename)
    cicj = cicjForInfiniteEnv(L);
    [I, Gplus, Gminus, Gp2, Gp3, GpGm, Gp2Gm, L2] = getRenyiNegNecessities(cicj, L/2 - l + 1, L/2, L/2 + l);
    alphas = -3.14*2:0.01:3.14*2;
    r1 = zeros(1, length(alphas));
    r2 = zeros(1, length(alphas));
    r3 = zeros(1, length(alphas));
    for i = 1:length(alphas)
        r1(i) = fluxRenyi1(I, Gplus, Gminus, alphas(i), L2);
        r2(i) = fluxRenyi2(I, Gplus, Gminus, alphas(i), L2);
        r3(i) = fluxRenyi3(I, Gplus, Gminus, Gp2, Gp3, GpGm, Gp2Gm, alphas(i), L2);
    end
    clear cicj
    save(filename)
end 