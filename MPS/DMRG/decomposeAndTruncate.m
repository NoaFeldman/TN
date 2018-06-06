function [psi, truncErr] = decomposeAndTruncate(M, k, psi, dir, opts)
    % perform an SVD decomposition on M.
    % Assign results to psi(k), psi(k+1). 
    [psi(k), psi(k+1), I] = myOrthoQS(M, [1, 2], dir, opts);
    truncErr = I.svd2tr;
    psi(k).info.itags(3) = strcat(int2str(k), 'a', psi(k).info.itags(3));
    psi(k+1).info.itags(1) = strcat(int2str(k), 'a', psi(k+1).info.itags(1));
end

function [l, r, I] = myOrthoQS(M, idx, dir, opts)
    if (length(M.Q) == 3)
        [r, l, I] = orthoQS(M, idx, dir, opts{:});
    else
        [l, r, I] = orthoQS(M, idx, dir, opts{:});
    end
    if (I.Nkeep == I.Ntot)
        return
    end
    originalNkeep = I.Nkeep;
    Nkeep = originalNkeep;
    while (Nkeep < I.Ntot && abs(I.svd(originalNkeep) - I.svd(Nkeep + 1)) / I.svd(originalNkeep) < 1e-2)
        Nkeep = Nkeep + 1;
    end
    if (Nkeep ~= originalNkeep)
        [l, r, I] = orthoQS(M, idx, dir, 'Nkeep', Nkeep);
    end
end

function [l, r] = sectorTrunc(M, dir, opts)
    r = QSpace;
    l = QSpace;
    Ms = getSectors(M);
    for i = 1:length(Ms)
        [currL, currR] = orthoQS(Ms(i), [1, 2], dir, opts{:});
        r = r + currR;
        l = l + currL;
    end
end

function Ms = getSectors(M)
    Ms = QSpace(length(M.data));
    for i = 1:length(Ms)
        Ms(i) = M;
        for j = 1:4
        Ms(i).Q{j} = M.Q{j}(i);
        end
        Ms(i).data = M.data(i);
    end
end
