function mirror = mirrorState(psi)
    % Assuming psi is left cannonicized.
    % Currently only working for U)1) symmetry.
    mirror = QSpace(length(psi));
    for i = 1:length(psi)
        mirrorInd = length(psi) + 1 - i;
        mirror(mirrorInd) = psi(i);
        mirror(mirrorInd).info.itags{3} = psi(i).info.itags{1};
        mirror(mirrorInd).info.itags{1} = psi(i).info.itags{3};
        mirror(mirrorInd).Q{3} = psi(i).Q{1};
        mirror(mirrorInd).Q{1} = psi(i).Q{3};
    end
end
    