function sNAFromAnal(L, tStep, tStepNum, filename)
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

    x = L / 4 - 2 : L / 4 + 2;
    s = zeros(length(x), tStepNum + 1);
    sFull = zeros(1, tStepNum + 1);
    U = realSpaceToDualSpace(L);
    for step = 0 : tStepNum
        [~, v] = eig(cicj(1:L/2, 1:L/2));
        for i = 1 : length(v)
            f(i) = v(i, i);
        end
%         res(step + 1)  = sum(getEE(f, x, L)); 
%         res(:, step + 1)  = getSNA(1, f, x, L);
        res(:, step + 1)  = getEE(f, x, L);
        sFull(step + 1) = sum(res(:, step + 1));
        ckcq = U' * cicj * U;
        ckcq = expectedCkCqMatrix(L, ckcq, tStep);
        cicj = U * ckcq * U';
    end 
    t = 0 : tStepNum;
    t = t * sFull;
    plot(t, res);
    legendInfo{1} = 'sFull';
    for i = 1 : length(x)
        plot(t(:), res(i, :));
        legendInfo{i + 1} = strcat('$s^z = $', num2str(L/4 - x(i)));
        hold on
    end
    legend(legendInfo, 'Interpreter', 'latex');
    set(gca, 'XScale', 'log');
    savefig(filename);
    hold off
end

function S = getEE(f, x, L)
    dn = 1e-3;
    snTop = getSNA(1 + dn / 2, f, x, L);
    snBottom = getSNA(1 - dn / 2, f, x, L);
    S = -(snTop - snBottom) / dn;
end