function M = getCiCjMatrix(psi, L)
    % Creates M_ij = <psi|Ci^dagger Cj|psi>
    % M will be an L*L matrix of the first L sites of the chain.
    M = zeros(length(psi), length(psi));
    for i = 1 : L
        for j = i : L
            M(i, j) = getCiCj(i, j, psi);
            M(j, i) = conj(M(i, j));
        end
    end
end