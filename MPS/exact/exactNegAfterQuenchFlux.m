function r3 = exactNegAfterQuenchFlux(L, u1, v1, v2, tFirstStep, tStep, tStepNum, directory, fileNameAddition, threadNum)
    maxNumCompThreads(threadNum);
    cicjA = getCiCj0Matrix(L/2);
    cicj = zeros(L);
    cicj(1:L/2, 1:L/2) = cicjA;
    cicj(L/2 + 1 : L, L/2 + 1 : L) = cicjA;
    alphas = -3.14:0.01:3.14;
    
    U = realSpaceToDualSpace(L);
    ckcq = U' * cicj * U;
    ckcq = expectedCkCqMatrix(L, ckcq, tFirstStep);
    cicj = U * ckcq * U';
    r1 = zeros(tStepNum + 1, length(alphas));
    r2 = zeros(tStepNum + 1, length(alphas));
    r3 = zeros(tStepNum + 1, length(alphas));
    for step = 0 : tStepNum
        [I, Gplus, Gminus, Gp2, Gp3, GpGm, Gp2Gm, L2] = getRenyiNegNecessities(cicj, u1, v1, v2);
        for i = 1:length(alphas)
            r1(step + 1, i) = fluxRenyi1(I, Gplus, Gminus, alphas(i), L2);
            r2(step + 1, i) = fluxRenyi2(I, Gplus, Gminus, alphas(i), L2);
            r3(step + 1, i) = fluxRenyi3(I, Gplus, Gminus, Gp2, Gp3, GpGm, Gp2Gm, alphas(i), L2);
        end
        ckcq = U' * cicj * U;
        ckcq = expectedCkCqMatrix(L, ckcq, tStep);
        cicj = U * ckcq * U';
    end 
    t = 0 : tStepNum;
    t = t * tStep + tFirstStep;
    save(strcat(directory, '/exactNegQuenchFlux', int2str(L), '_', int2str(v1 - u1 + 1), '-', int2str(v2 - v1), '_', fileNameAddition)...
        , 't', 'r1', 'r2', 'r3', 'u1', 'v1', 'v2', 'alphas');
end
