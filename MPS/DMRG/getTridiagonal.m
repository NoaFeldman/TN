function [T, base] = getTridiagonal(HR, HL, HR2, HL2, H, k, psi)
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
    accuarcy = 1e-10; % 1e-12;
    if length(psi(1).info.itags) == 3
        bc = 'open';
    else
        bc = 'periodic';
    end
    
    if strcmp(bc, 'open')
        v = contract(psi(k), 3, psi(k+1), 1);
        % Small innaccuracies ruin everything!
        v = v / sqrt(getscalar(contract(v, '1234', v, '1234*')));
        base(1) = v;
        Hv = applyHToM(HR, HL, HR2, HL2, H, v, k, bc);
        alpha = getscalar(contract(v, '1234', Hv, '1234*'));
        w = Hv - alpha * v;
        beta = sqrt(getscalar(contract(w, '1234', w, '1234*')));
    else
        v = contract(psi(k), 4, psi(k+1), 1);
        % Small innaccuracies ruin everything!
        base(1) = v;
        Hv = applyHToM(HR, HL, HR2, HL2, H, v, k, bc);
        alpha = getscalar(contract(v, '123456', Hv, '123456*'));
        w = Hv - alpha * v;
        beta = sqrt(getscalar(contract(w, '123456', w, '123456*')));
    end
    T(1, 1) = alpha;
    counter = 1;
    formBeta = 2*beta; % This is just some value to init formBeta > beta.
    while beta > accuarcy & counter <= 50 & beta < formBeta
        T(counter, counter+1) = beta;
        T(counter+1, counter) = beta;
        counter = counter + 1;
        v = w / beta;
        base(counter) = v;
        Hv = applyHToM(HR, HL, HR2, HL2, H, v, k, bc);
        if strcmp(bc, 'open')
            alpha = getscalar(contract(v, '1234', Hv, '1234*'));
        else
            alpha = getscalar(contract(v, '123456', Hv, '123456*'));            
        end
        T(counter, counter) = alpha;
        w = Hv - alpha * v - beta * base(counter - 1);
        formBeta = beta;
        if strcmp(bc, 'open')
            beta = sqrt(getscalar(contract(w, '1234', w, '1234*')));
        else
            beta = sqrt(getscalar(contract(w, '123456', w, '123456*')));
        end
    end
end
    

        