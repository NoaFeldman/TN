function cnPsi = applyCn(psi, n)
    % (relevant for non interacting spin 1/2 systems)
    % Assuming psi is left canonical.
    % Applys cn = exp(i * pi * sum_k=1:n-1 { (s^z_k+1) / 2 }) * s^-_n
    cnPsi = QSpace(length(psi));
    [S,IS]=getLocalSpace('Spin',0.5,'-A');
    % exp(i * pi * (s^z_i + 1) / 2) = - s^z_i
    signOp = -2 * S(1);
    signOp.Q = signOp.Q(1:2);
    signOp.info.itags = signOp.info.itags(1:2);
    for i = 1 : n - 1
        cnPsi(i) = contract(signOp, 2, psi(i), 2, [2 1 3]);
        cnPsi(i).info.itags{2} = strcat(int2str(i), 's');
    end
    s = sqrt(2) * S(2);
    s.info.itags = {strcat(int2str(n), 's'), strcat(int2str(n), 's*'), 'm*'};
    cnPsi(n) =  contract(s, 2, psi(n), 2, [3 1 4 2]);
    for i = n + 1 : length(psi)
        cnPsi(i) = psi(i);
    end
end