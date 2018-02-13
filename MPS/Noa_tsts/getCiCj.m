function cicj = getCiCj(i, j , psi) 
    % Assuming j >= i
    CStates = QSpace(length(psi), length(psi));
    for n = 1 : length(CStates)
        CStates(n, :) = applyCn(psi, n);
    end
    res = QSpace;
    CiStates = CStates(i, :);
    CjStates = CStates(j, :);
    if (i < j)
        res = contract(CjStates(i), '12', CiStates(i), '12*');
        for n = i+1 : j-1
            res = contract(contract(res, 1, CjStates(n), 1), '13', ...
                CiStates(n), '12*', [2 3 1]);
        end
        res = contract(contract(res, '13', CjStates(j), '14'), '12', ...
            CiStates(j), '12*');
    else % i <= j according to getCiCjMatrix
        res = contract(CjStates(i), '124', CiStates(i), '124*');
    end
    for n = j+1 : length(CiStates)
        res = contract(contract(res, 1, CjStates(n), 1), '12', ...
            CiStates(n), '12*');
    end
    cicj = getscalar(contract(res, '12', getIdentity(res, 1), '21'));
end         