function [psi, truncErr] = decomposeAndTruncate(M, k, psi, dir, opts)
    % perform an SVD decomposition on M.
    % Assign results to psi(k), psi(k+1).    
    [psi(k), psi(k+1), I] = orthoQS(M, [1, 2], dir, opts{:});
    disp(strcat('k = ', int2str(k)));
    disp('M = ');
    disp(M);
    disp('psi(k) = ');
    disp(psi(k));
    disp('psi(k+1) = ');
    disp(psi(k+1));
    truncErr = I.svd2tr;
    psi(k).info.itags(3) = strcat(int2str(k), 'a', psi(k).info.itags(3));
    psi(k+1).info.itags(1) = strcat(int2str(k), 'a', psi(k+1).info.itags(1));
end