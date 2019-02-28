function [psi, truncErr] = decomposeAndTruncate(M, k, psi, dir, opts)
    % perform an SVD decomposition on M.
    % Assign results to psi(k), psi(k+1).
    if length(psi(1).info.itags) == 3
        bc = 'open';
    else
        bc = 'periodic';
    end
    if strcmp(bc, 'open')
        [psi(k), psi(k+1), I] = myOrthoQS(M, [1, 2], dir, opts);
    else
        [psi(k), psi(k+1), I] = myOrthoQS(M, 1:3, dir, opts);
    end
    truncErr = I.svd2tr;
    psi(k).info.itags(length(psi(k).info.itags)) = ...
        strcat(int2str(k), 'a', psi(k).info.itags(length(psi(k).info.itags)));
    psi(k+1).info.itags(1) = strcat(int2str(k), 'a', psi(k+1).info.itags(1));
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
