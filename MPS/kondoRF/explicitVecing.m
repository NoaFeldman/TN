function [LiouMat, rhoVec] = explicitVecing(Liou, rho)
    % Following https://arxiv.org/pdf/1510.08634.pdf
    % I don't use this method anymore and define the Liouvillian as amatrix
    % using the TN formalism, but this is still useful for debugging
    % reasons.
    HMat = blkdiag(Liou.H);
    rhoMat = blkdiag(rho);
    n = length(rhoMat);
    rhoVec = reshape(rhoMat, [n^2, 1]);
    Liou.L = 1e-99*Liou.H + Liou.L;
    Liou.LdagL = 1e-99*Liou.H + Liou.LdagL;
    % Modify this for left diagonalization!
    LdagLMat = blkdiag(Liou.LdagL);
    LMat = blkdiag(Liou.L);
        
    LiouMat = -1i * (kron(eye(n), HMat) - kron(transpose(HMat), eye(n))) + ...
              kron(conj(LMat), LMat) - ...
              0.5 .* (kron(eye(n), LdagLMat) + kron(transpose(LdagLMat), eye(n)));
end