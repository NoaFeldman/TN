function folded = foldChain(psi)
    % folds      _
    % |>- -|>- -|_|- -<|
    % to  _     _
    % |>-| |   | |-|>__
    %    | |- -| |  _  |
    % |>-|_|   |_|-|_|_|
    folded = QSpace(length(psi) / 2);
    k = length(psi);
    while k > length(psi) / 2
        [psi, k] = shiftWorkingSite(psi, k, '<<');
    end
    for lowIndex = k : -1 : 1
        highIndex = length(psi) + 1 - lowIndex;
        sID = getIdentity(real(psi(lowIndex)), 2, real(psi(highIndex)), 2, 'o');
        sID.info.itags{3} = strcat('s', int2str(lowIndex), '-', int2str(highIndex), '*');
        aLeftID = getIdentity(real(psi(lowIndex)), 1, real(psi(highIndex)), 3, 'o');
        aLeftID.info.itags{3} = strcat('a', int2str(lowIndex - 1), '-', int2str(highIndex), '*');
        if lowIndex < k
            aRightID = getIdentity(real(psi(lowIndex)), 3, real(psi(highIndex)), 1, 'o');
            aRightID.info.itags{3} = strcat('a', int2str(lowIndex), '-', int2str(highIndex - 1), '*');
            folded(lowIndex) = contract(...
                contract(psi(lowIndex), 3, aRightID, 1), 3, ...
                psi(highIndex), 1, [1 5 2 4 3]);
        else
            folded(lowIndex) = contract(psi(lowIndex), 3, psi(highIndex), 1, [1 4 2 3]);
        end
        folded(lowIndex) = contract(aLeftID, '12*', folded(lowIndex), '12');
        if lowIndex < k
            folded(lowIndex) = contract(folded(lowIndex), '23', sID, '12*', [1 3 2]);
        else
            folded(lowIndex) = contract(folded(lowIndex), '23', sID, '12*');
        end
    end  
end
    
    