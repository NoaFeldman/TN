function n = getPsiNorm(psi, k)
    % k is the index of the site in mixed canonical form.
    nLeft = getIdentity(psi(1), 1, [2 1]);
    for i = 1:k
        nLeft = contract(contract(nLeft, 1, psi(i), 1), '12', psi(i), '12*');
    end
    nRight = getIdentity(psi(length(psi)), 3, [2 1]);
    for i = length(psi):-1:k+1
        nRight =  contract(contract(nRight, 1, psi(i), 3), '13', psi(i), '32*');
    end
    n = getscalar(contract(nLeft, '12', nRight, '12'));
end