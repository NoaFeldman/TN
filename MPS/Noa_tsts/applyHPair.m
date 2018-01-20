function psi = applyHPair(trotterGates, k, psi, dir, orthoOpts)
    % apply the HPairs(k) to psi(k)---psi(k+1).
    % Apply SVD on the result and assign to psi(k) and psi(k+1).
    M = contract(psi(k), 3, psi(k+1), 1);
    % trotterGate itags example : { 3s*, 4s*, 3s, 4s }
    % M itag example: { 3s, 4s, 34s* }
    M = contract(M, '23', trotterGates(k), '12', [1 3 4 2]);
    psi = decomposeAndTruncate(M, k, psi, dir, orthoOpts);
end
    