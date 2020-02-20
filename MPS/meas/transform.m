function transformed = transform(psi, copyNum)
    copies = QSpace(copyNum, length(psi));
    omega = exp(1i * 2*pi / copyNum);
    for j = 1:length(psi)
        emptyCopies = getEmptyCopies(copies, j, copyNum);
        for combination = 0 : 2^copyNum-1
            originalCopies = getCopiesInCombination(psi(j), copyNum, combination);
            
        end
    end
    
    transformed = psi;
end

function emptyCopies = getEmptyCopies(copies, siteNum, copyNum)
    [S,IS]=getLocalSpace('Spin',0.5,'-A');
    emptyCopies = QSpace(copyNum);
    if siteNum = 1
        for n = 1:copyNum
            epmtyCopies(n) = getIdentity(S(1), 3, IS.E, 1);
        end
    else
        for n = 1:copyNum
            epmtyCopies(n) = getIdentity(copies(n, siteNum - 1), 3, IS.E, 1);
            if siteNum == length(copies(1, :))
                epmtyCopies(n).Q{3} = -epmtyCopies(n).Q{3};
                epmtyCopies(n).info.itags{3} = strcat(int2str(siteNum), 'a');
            end
        end
    end
end

function copies = getCopiesInCombination(site, copyNum, combination)
    copies = QSpace(copyNum);
    [S,IS]=getLocalSpace('Spin',0.5,'-A');
    IS.E.info.itags = {site.info.itags{2}, strcat(site.info.itags{2}, '*')};
    projUp = IS.E;
    projUp.Q{1} = projUp.Q{1}(2);
    projUp.Q{2} = projUp.Q{2}(2);
    projUp.data = projUp.data(2);
    projDown = IS.E;
    projDown.Q{1} = projDown.Q{1}(1);
    projDown.Q{2} = projDown.Q{2}(1);
    projDown.data = projDown.data(1);
    for n = 1:copyNum
        population = (double(bitand(combination, 2^(n-1)) > 0) - 0.5)*2;
        if population == 1
            projector = projUp;
        else
            projector = projDown;
        end
        copies(n) = contract(site, 2, projector, 2, [1 3 2]);
    end
end