function [psi, truncErr] = decomposeAndTruncate(M, k, psi, dir, opts)
    % perform an SVD decomposition on M.
    % Assign results to psi(k), psi(k+1). 
    M = flip(M);
    [psi(k), psi(k+1), I] = orthoQS(M, [1, 2], dir, opts{:});
    truncErr = I.svd2tr;
    psi(k).info.itags(3) = strcat(int2str(k), 'a', psi(k).info.itags(3));
    psi(k+1).info.itags(1) = strcat(int2str(k), 'a', psi(k+1).info.itags(1));
end

function flipped = flip(M)
    flipped = M;
    for i = 1 : length(M.data)/2 + 1
        j = length(M.data) + 1 - i;
        flipped.data{i} = M.data{j};
        flipped.data{j} = M.data{i};
        for k = 1: length(M.Q)
            flipped.Q{k}(i) = M.Q{k}(j);
            flipped.Q{k}(j) = M.Q{k}(i);
        end
    end
end
    