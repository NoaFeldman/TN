function psi = coupleStates(psiL, psiR)
    % Create a state of the two chains described by psiL and psiR
    % Assuming both states psiL and psiR are left canonitized.
    psi = QSpace(length(psiL) + length(psiR));
    psi(1:length(psiL)) = psiL;
    psi(length(psiL)) = contract(psi(length(psiL)), 3, ...
        getIdentity(psi(length(psiL)), 3, '-0'), '1*');
    for i = 1 : length(psiR)
        index = length(psiL) + i;
        psi(index) = psiR(i);
        psi(index).info.itags = {strcat(int2str(index - 1), 'a'), ...
                                 strcat(int2str(index), 's'), ...
                                 strcat(int2str(index), 'a*')};
        if (index == length(psi))
            psi(index).info.itags{3} = strcat(int2str(index), 'a');
        end
    end
end