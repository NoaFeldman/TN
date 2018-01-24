function [res, t] = getCK1CK2(psi, nk1, nk2)
% For non interacting spin chains.
% TODO
    tic;
    k1 = nk1 * 2 * pi / (length(psi) + 1);
    k2 = nk2 * 2 * pi / (length(psi) + 1);
    res = 0;
    cn = QSpace(length(psi), length(psi));
    for n = 1 : length(psi)
        cn(n, :) = applyCn(psi, n);
        t(n) = toc;
    end    
    for n = 1 : length(psi)
        for m = 1 : length(psi)
            res = res + cos(k1 * n) * cos(k2 * m) * ...
                getOverlap(cn(n, : ), cn(m, :), length(psi));
        end
        t(n + length(psi)) = toc;
    end
end