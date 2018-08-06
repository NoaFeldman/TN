function sNAExact(L, u, v, tFirstStep, tStep, tStepNum, fileNameAddition, alphas)
    path(path, [pwd, '/theo']);
    
    cicjA = getCiCj0Matrix(L/2);
    cicj = zeros(L, L);
    cicj(1:L/2, 1:L/2) = cicjA;
    cicj(L/2 + 1 : L, L/2 + 1 : L) = cicjA;
    
    x = (v-u) / 2 - 10 : 2 : (v-u) / 2 + 10;
    s = zeros(length(x), tStepNum + 1);
    s1 = zeros(length(x), tStepNum + 1);
    s2 = zeros(length(x), tStepNum + 1);
    s3 = zeros(length(x), tStepNum + 1);
    s4 = zeros(length(x), tStepNum + 1);
    s5 = zeros(length(x), tStepNum + 1);
    s1Alpha = zeros(length(alphas), tStepNum + 1);
    s2Alpha = zeros(length(alphas), tStepNum + 1);
    s3Alpha = zeros(length(alphas), tStepNum + 1);
    s4Alpha = zeros(length(alphas), tStepNum + 1);
    s5Alpha = zeros(length(alphas), tStepNum + 1);
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
        s1(:, step + 1)  = getSNA(1, f, x, L);
        s2(:, step + 1)  = getSNA(2, f, x, L);
        s3(:, step + 1)  = getSNA(3, f, x, L);
        s4(:, step + 1)  = getSNA(4, f, x, L);
        s5(:, step + 1)  = getSNA(5, f, x, L);
        s1Alpha(:, step + 1) = getSAlpha(1, alphas, f, 1);
        s2Alpha(:, step + 1) = getSAlpha(2, alphas, f, 1);
        s3Alpha(:, step + 1) = getSAlpha(3, alphas, f, 1);
        s4Alpha(:, step + 1) = getSAlpha(4, alphas, f, 1);
        s5Alpha(:, step + 1) = getSAlpha(5, alphas, f, 1);
        s(:, step + 1)  = getEE(f, x, L);
        sFull(step + 1) = sum(s(:, step + 1));
        
        ckcq = U' * cicj * U;
        ckcq = expectedCkCqMatrix(L, ckcq, tStep);
        cicj = U * ckcq * U';
    end 
    t = 0 : tStepNum;
    t = t * tStep;
    save(strcat('theoSP', int2str(L), '_', int2str(u), '-', int2str(v), '_', num2str(tStep*tStepNum), '_', fileNameAddition)...
        , 't', 's', 's1', 's2', 's3', 's4', 's5', ...
        'alphas', 's1Alpha', 's2Alpha', 's3Alpha', 's4Alpha', 's5Alpha', 'sFull');
end
