L = 10000;
lRights = 200:50:500;
lLeft = 200;
cicj = getCiCj0Matrix(L);
alphas = -3.14:1e-2:3.14;
negs = zeros(length(ls), length(alphas));
for i = 1:length(lRights)
    lRight = lRights(i);
    [I, Gplus, Gminus, Gp2, Gp3, GpGm, Gp2Gm, L2] = getRenyiNegNecessities(cicj, L/2 - lLeft  +1, L/2, L/2 + lRight);
    for j = 1:length(alphas)
        alpha = alphas(j);
        negs1(i, j) = fluxRenyi1(I, Gplus, Gminus, alpha, L2);
        negs2(i, j) = fluxRenyi2(I, Gplus, Gminus, alpha, L2);
        negs3(i, j) = fluxRenyi3(I, Gplus, Gminus, alpha, L2);
    end
end
save('negtst');
    