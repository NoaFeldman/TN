function [HL, HR, psi, E0, k, M] = dmrgStep(HL, HR, H, psi, k, dir, opts)
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
    [M, E0] = lanczos(HL(k1), HR(k2), H, k1, psi);
%     [psi(k1), psi(k2)] = orthoQS(M, [1, 2], dir, 'Nkeep', Nkeep);
%     psi(k1).info.itags(3) = strcat(int2str(k1), 'a', psi(k1).info.itags(3));
%     psi(k2).info.itags(1) = strcat(int2str(k1), 'a', psi(k2).info.itags(1));
    psi = decomposeAndTruncate(M, k1, psi, dir, opts);
    if (strcmp(dir, '>>'))
        HL(k+1) =  getHLR(H, psi, k, '>>', HL(k));
        k = k+1;
    else
        HR(k-1) =  getHLR(H, psi, k, '<<', HR(k));
        k = k - 1;
    end
end
