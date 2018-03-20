function sNAFromAnal(L, tStep, tStepNum, filename)
    cicj = getCiCj0(L);
    x = L / 4 - 5 : L / 4 + 5;
    res = zeros(length(x), tStepNum + 1);
    for step = 0 : tStepNum
        [~, v] = eig(cicj(1:L/2, 1:L/2));
        for i = 1 : length(v)
            f(i) = v(i, i);
        end
        res(:, step + 1)  = getSNA(1, f, x);
        S = realSpaceToDualSpace(L);
        ckcq = S' * cicj * S;
        ckcq = expectedCkCqMatrix(L, ckcq, tStep);
        cicj = S * ckcq * S';
    end 
    t = 0 : tStepNum;
    t = t * tStep;
    for i = 1 : length(x)
        scatter(t(:), res(i, :));
        legendInfo{i} = strcat('$N_A = $', num2str(x(i)));
        hold on
    end
    legend(legendInfo, 'Interpreter', 'latex');
    savefig(filename);
    hold off
end

function C = getCiCj0(L)
C = zeros(L, L);
    for i = 1 : L/2
        for j = 1 : L/2
            if (i == j) 
                C(i, j) = 0.5;
            else
                C(i, j) = sin(pi * (i-j) / 2) / (pi * (i - j));
                C(j, i) = C(i, j);
            end
        end
    end
    for i = L/2 + 1 : L
        for j = L/2 + 1 : L
            if (i == j) 
                C(i, j) = 0.5;
            else
                C(i, j) = sin(pi * (i-j) / 2) / (pi * (i - j));
                C(j, i) = C(i, j);
            end
        end
    end
end