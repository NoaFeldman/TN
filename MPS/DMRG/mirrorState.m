function mirror = mirrorState(psi)
    % Assuming psi is left cannonicized.
    % Not handling mirror.info.itags.
    % Currently only works 
    mirror = QSpace(length(psi));
    for i = 1:length(psi)
        mirrorInd = length(psi) + 1 - i;
        mirror(mirrorInd) = permute(psi(i), length(psi(i).Q):-1:1);
    end
end
    