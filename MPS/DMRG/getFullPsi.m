function p = getFullPsi(psi)
    p = psi(1);
    for i = 2 : length(psi)
        p = contract(p, i + 1, psi(i), 1);
    end
end