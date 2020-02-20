function mergedCopies = toCopies(psi, n)
    copies = QSpace(n, length(psi));
    for i = 1:n
        copies(i, :) = psi;
        for j = 1:length(psi)
            copies(i, j) = permute(copies(i, j), [1 3 2]);
            copies(i, j).info.itags{3} = strcat(copies(i, j).info.itags{3}, 'c', int2str(i));
        end
    end
    for j = 1:length(psi)
        mergedCopies(j) = merge(copies(:, j));
    end
end

function merged = merge(copies)
    idIn = getIdentity(copies(1), 1, copies(2), 1, copies(1).info.itags{1});
    idOut = getIdentity(copies(1), 2, copies(), 2, copies(2).info.itags{2});
    merged = contract(contract(idIn, '1*', copies(1), 1), 1, copies(2), 1);
    if ~isempty(strfind(merged.info.itags{2}, '*')) % Mid chain (for OC on the right)
        merged = contract(merged, '24', ...
            idOut, '12');
    else
        merged = contract(merged, '24', ...
            idOut, '12*');
    end
end        