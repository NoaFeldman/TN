function psi = applyHPair(HPairs, k, psi, dir, orthoOpts)
    % apply the HPairs(k) to psi(k)---psi(k+1).
    % Apply SVD on the result and assign to psi(k) and psi(k+1).
    M = contract(psi(k), 3, psi(k+1), 1);
    % HPairs itags example: { 2s*, 2s, 3s*, 3s }
    % M itags example: { 1a*, 2s, 3s, 3a }
    M = contract(M, '23', HPairs(k), '13', [1 3 4 2]);
    psi = decomposeAndTruncate(M, k, psi, dir, orthoOpts);
end
    