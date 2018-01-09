function [T, base] = getTridiagonal(HL, HR, H, k, psi)
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
    accuarcy = 0.001;
    v = contract(psi(k), 3, psi(k+1), 1);
    % Small innaccuracies ruin everything!
    v = v / sqrt(getscalar(contract(v, '1234', v, '1234*')));
    base(1) = v;
    Hv = applyHToM(HR, HL, H, v, k);
    alpha = getscalar(contract(v, '1234', Hv, '1234*'));
    T(1, 1) = alpha;
    w = Hv - alpha * v;
    beta = sqrt(getscalar(contract(w, '1234', w, '1234*')));
    counter = 1;
    while beta > accuarcy 
        disp(['k = ', num2str(k) , ', beta = ', num2str(beta)]);
        T(counter, counter+1) = beta;
        T(counter+1, counter) = beta;
        counter = counter + 1;
        v = w / beta;
        base(counter) = v;
        Hv = applyHToM(HR, HL, H, v, k);
        alpha = getscalar(contract(v, '1234', Hv, '1234*'));
        T(counter, counter) = alpha;
        w = Hv - alpha * v - beta * base(counter - 1);
        beta = sqrt(getscalar(contract(w, '1234', w, '1234*')));
    end

        