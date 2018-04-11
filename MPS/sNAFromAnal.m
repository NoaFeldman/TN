function sNAFromAnal(L, tStep, tStepNum)
    path(path, [pwd, '/MPSPACK_v3.0']);
    path(path, [pwd, '/Noa_tsts']);
    startup;
    tic;
    % TODO take to outside func
    ckcqA = zeros(L/2, L/2);
    for i = 1 : L/4 % L/4 + 1 : L/2
        ckcqA(i, i) = 1;
    end
    SA = realSpaceToDualSpace(L/2); 
    cicjA = SA * ckcqA * SA';
    cicj = zeros(L, L);
    cicj(1:L/2, 1:L/2) = cicjA;
    cicj(L/2 + 1 : L, L/2 + 1 : L) = cicjA;
    disp('got cicj');
    toc;
    
    x = L / 4 - 2 : L / 4 + 2;
    s = zeros(length(x), tStepNum + 1);
    p = zeros(length(x), tStepNum + 1);
    sFull = zeros(1, tStepNum + 1);
    U = realSpaceToDualSpace(L);
    disp('got U');
    toc;
    for step = 0 : tStepNum
        [~, v] = eig(cicj(1:L/2, 1:L/2));
        for i = 1 : length(v)
            f(i) = v(i, i);
        end
        disp(strcat('step = ', num2str(step), ', line 31'));
        toc;
        p(:, step + 1)  = getSNA(1, f, x, L);
        disp(strcat('step = ', num2str(step), ', line 34'));
        toc;
        s(:, step + 1)  = getEE(f, x, L);
        sFull(step + 1) = sum(s(:, step + 1));
        disp(strcat('step = ', num2str(step), ', line 38'));
        toc;
        ckcq = U' * cicj * U;
        ckcq = expectedCkCqMatrix(L, ckcq, tStep);
        cicj = U * ckcq * U';
        disp(strcat('step = ', num2str(step)));
        toc;
        if (mod(step, 100) == 0)
                save(strcat('theoSP', int2str(L), 'step_', num2str(step)), 's', 'p', 'sFull', 't');
        end
    end 
    t = 0 : tStepNum;
    t = t * tStep;
    save(strcat('theoSP', int2str(L), '_', num2str(tStepNum)), 's', 'p', 'sFull', 't');
end

function S = getEE(f, x, L)
    dn = 1e-3;
    snTop = getSNA(1 + dn / 2, f, x, L);
    snBottom = getSNA(1 - dn / 2, f, x, L);
    S = -(snTop - snBottom) / dn;
end