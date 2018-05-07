function psi = coupleStates(psiL, psiR)
    % Create a state of the two chains described by psiL and psiR
    % Assuming both states psiL and psiR are left canonitized.
    psi = QSpace(length(psiL) + length(psiR));
    psi(1:length(psiL)) = psiL;
    psi(length(psiL)) = contract(psi(length(psiL)), 3, ...
        getIdentity(psi(length(psiL)), 3, '-0'), '1*');
    psiR = mirrorState(psiR);
    for i = 1 : length(psiR)
        index = length(psiL) + i;
        psi(index) = psiR(i);
        psi(index).info.itags = {strcat(int2str(index - 1), 'a*'), ...
                                 strcat(int2str(index), 's'), ...
                                 strcat(int2str(index), 'a')};
        if (i == 1)
            psi(index).info.itags{1} = strcat(int2str(index-1), 'a');
        end
    end
    workingSiteIndex = length(psiL) + 1;
    while workingSiteIndex < length(psi)
        [psi, workingSiteIndex] = shiftWorkingSite(psi, workingSiteIndex, '>>');
    end
end
