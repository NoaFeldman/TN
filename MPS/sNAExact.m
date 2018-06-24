function sNAExact(L, LA, tFirstStep, tStep, tStepNum, fileNameAddition)
    path(path, [pwd, '/theo']);
    
    cicjA = getCiCj0Matrix(L/2);
    cicj = zeros(L, L);
    cicj(1:L/2, 1:L/2) = cicjA;
    cicj(L/2 + 1 : L, L/2 + 1 : L) = cicjA;
    
    x = LA / 2 - 10 : 2 : LA / 2 + 10;
    s = zeros(length(x), tStepNum + 1);
    p = zeros(length(x), tStepNum + 1);
    sFull = zeros(1, tStepNum + 1);
    
    % Sanity check - S_A == S_B?
    sb = zeros(length(x), tStepNum + 1);
    pb = zeros(length(x), tStepNum + 1);
    sFullb = zeros(1, tStepNum + 1);
    
    U = realSpaceToDualSpace(L);
    ckcq = U' * cicj * U;
    ckcq = expectedCkCqMatrix(L, ckcq, tFirstStep);
    cicj = U * ckcq * U';
    for step = 0 : tStepNum
        [~, v] = eig(cicj(1:LA, 1:LA));
        for i = 1 : length(v)
            f(i) = v(i, i);
        end
        p(:, step + 1)  = getSNA(1, f, x, L);
        s(:, step + 1)  = getEE(f, x, L);
        sFull(step + 1) = sum(s(:, step + 1));
        
        [~, v] = eig(cicj(LA+1:L, LA+1:L));
        for i = 1 : length(v)
            f(i) = v(i, i);
        end
        pb(:, step + 1)  = getSNA(1, f, L/2 - x, L);
        sb(:, step + 1)  = getEE(f, L/2 - x, L);
        sFullb(step + 1) = sum(s(:, step + 1));
        
        
        ckcq = U' * cicj * U;
        ckcq = expectedCkCqMatrix(L, ckcq, tStep);
        cicj = U * ckcq * U';
%         if (mod(step, 100) == 0)
%             t = 0 : step;
%             t = t * tStep;
%             save(strcat('theoSP', int2str(L), '_', int2str(LA), '_step_', num2str(step)), 's', 'p', 'sFull', 't');
%         end
    end 
    t = 0 : tStepNum;
    t = t * tStep + tFirstStep;
    save(strcat('theoSP', int2str(L), '_', int2str(LA), '_', num2str(tStepNum), '_', fileNameAddition)...
        , 't', 's', 'p', 'sFull', 'sb', 'pb', 'sFullb');
end

