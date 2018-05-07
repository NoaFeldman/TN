function cicj = getCiCj(i, j , psi) 
    % Assuming j >= i
    res = QSpace;
    CiStates = applyCn(psi, i);
    CjStates = applyCn(psi, j);
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
    if (isempty(res))
        cicj = 0;
    else
        cicj = getscalar(contract(res, '12', getIdentity(res, 1), '21'));
    end
end         