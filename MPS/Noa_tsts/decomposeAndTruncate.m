function psi = decomposeAndTruncate(M, k, psi, dir, orthoOpts)
    % perform an SVD decomposition on M.
    % Assign results to psi(k), psi(k+1).    
    [psi(k), psi(k+1)] = orthoQS(M, [1, 2], dir, orthoOpts{:});
    psi(k).info.itags(3) = strcat(int2str(k), 'a', psi(k).info.itags(3));
    psi(k+1).info.itags(1) = strcat(int2str(k), 'a', psi(k+1).info.itags(1));
end