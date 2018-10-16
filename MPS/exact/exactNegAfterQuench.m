function exactNegAfterQuench(L, u1, v1, v2, tFirstStep, tStep, tStepNum, directory, fileNameAddition)
    cicjA = getCiCj0Matrix(L/2);
    cicj = zeros(L);
    cicj(1:L/2, 1:L/2) = cicjA;
    cicj(L/2 + 1 : L, L/2 + 1 : L) = cicjA;
    Qs = -12:12;
    
    U = realSpaceToDualSpace(L);
    ckcq = U' * cicj * U;
    ckcq = expectedCkCqMatrix(L, ckcq, tFirstStep);
    cicj = U * ckcq * U';
    r3 = zeros(tStepNum + 1, length(Qs));
    for step = 0 : tStepNum
        for i = 1:length(Qs)
            r3(step + 1, i) = chargeRenyi3(cicj, Qs(i), u1, v1, v2);
        end
        ckcq = U' * cicj * U;
        ckcq = expectedCkCqMatrix(L, ckcq, tStep);
        cicj = U * ckcq * U';
    end 
    t = 0 : tStepNum;
    t = t * tStep + tFirstStep;
    save(strcat(directory, '/exactNegQuench', int2str(L), '_', int2str(v1 - u1 + 1), '-', int2str(v2 - v1), '_', fileNameAddition)...
        , 't', 'r3', 'u1', 'v1', 'v2', 'Qs');
end