function T = getTridiagonal(HL, HR, H, k, psi, dir)
    % Uses the Lanczos algorithm to compute a tridiagonal matrix based on
    % H.
    % HL, HR represent the two blocks of <H> from both sides, H(k) and
    % H(k+dir) the MPO of the kth and (k+dir)th sites, and psi(k) + psi(k+dir)
    % are the matrices to be optimized.
    % If dir == '>>' then the second MPS is on the left (i.e. k, k+1),
    % otherwise we use k, k-1. TODO
    %  __                       __
    % |  |---               ---|  |
    % |HL|    _|__   __|___    |HR|
    % |  |---|H(k)|-|H(k+1)|---|  |
    % |  |     |       |       |  |
    % |  |                     |  |
    % |  |  ___|___  __|_____  |  |
    % |__|-|psi(k)|-|psi(k+1)|-|__|
    % T is returned as an array of {alpha1, beta1; alpha2, beta2; ...}
    
    n = 4;
    T = double(n, 2);
    
    vcurr = contract(psi(k), 3, psi(k+1), 1);
    vprev = QSpace;
    beta = 1;
    for i = 1 : n
        w = apllyH(HR, HL, H, k, vcurr);
        alpha = contract(w, '1234*', vcurr, '1234').data{1};
        w = w - alpha * vcurr - beta*vprev;
        beta = contract(w, '1234*', w, '1234').data{1};
        T(i, 1) = alpha;
        T(i, 2) = beta;
        vprev = vcurr;
        vcurr = w / beta;
    end

    function w = apllyH(HR, HL, H, k, v)
        