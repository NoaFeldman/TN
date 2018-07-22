function sNAExact(L, u, v, tFirstStep, tStep, tStepNum, fileNameAddition)
    path(path, [pwd, '/theo']);
    
    cicjA = getCiCj0Matrix(L/2);
    cicj = zeros(L, L);
    cicj(1:L/2, 1:L/2) = cicjA;
    cicj(L/2 + 1 : L, L/2 + 1 : L) = cicjA;
    
    x = (v-u) / 2 - 10 : 2 : (v-u) / 2 + 10;
    s = zeros(length(x), tStepNum + 1);
    p = zeros(length(x), tStepNum + 1);
    sFull = zeros(1, tStepNum + 1);
    
  
    U = realSpaceToDualSpace(L);
    ckcq = U' * cicj * U;
    ckcq = expectedCkCqMatrix(L, ckcq, tFirstStep);
    cicj = U * ckcq * U';
    for step = 0 : tStepNum
        [~, V] = eig(cicj(u+1:v, u+1:v));
        for i = 1 : length(V)
            f(i) = V(i, i);
        end
        p(:, step + 1)  = getSNA(1, f, x, L);
        s(:, step + 1)  = getEE(f, x, L);
        sFull(step + 1) = sum(s(:, step + 1));
        
        ckcq = U' * cicj * U;
        ckcq = expectedCkCqMatrix(L, ckcq, tStep);
        cicj = U * ckcq * U';
    end 
    t = 0 : tStepNum;
    t = t * tStep;
    save(strcat('theoSP', int2str(L), '_', int2str(u), '-', int2str(v), '_', num2str(tStep*tStepNum), '_', fileNameAddition)...
        , 't', 's', 'p', 'sFull');
end
