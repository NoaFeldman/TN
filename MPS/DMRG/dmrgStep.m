function [HL, HR, HL2, HR2, psi, E0, k, M] = dmrgStep(HL, HR, HL2, HR2, H, psi, k, dir, opts)
    % Perform a single DMRG step:
    % 1. Contracts psi(k) and psi(k + dir) to get M.
    % 2. Performs lancsoz and get a new contracted M.
    % 3. Performs an SVD in order to get the new working site, at k + dir.
    % 4. Calculates HL(k) / HR(k) (according do dir)
    if (strcmp(dir, '>>'))
        k1 = k;
        k2 = k+1;
    else
        k1 = k -1;
        k2 = k;
    end
    [M, E0] = lanczos(HR(k2), HL(k1), HR2(k2), HL2(k1), H, k1, psi);
    [psi, truncErr] = decomposeAndTruncate(M, k1, psi, dir, opts);
    if (strcmp(dir, '>>'))
        [HL(k+1), HL2(k+1)] =  getHLR(H, psi, k, '>>', HL(k), HL2(k));
        k = k+1;
    else
        [HR(k-1), HR2(k-1)] =  getHLR(H, psi, k, '<<', HR(k), HR2(k));
        k = k - 1;
    end
end
