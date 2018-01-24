function cnPsi = applyCn(psi, n)
    % (relevant for non interacting spin 1/2 systems)
    % Assuming psi is left canonical.
    % Applys cn = exp(i * pi * sum_k=1^n-1 (s+_n s-_n)) s-_n
    cnPsi = QSpace(length(psi));
    [S,IS]=getLocalSpace('Spin',0.5,'-A');
    signOp = -2 * S(1);
    signOp.Q = signOp.Q(1:2);
    signOp.info.itags = signOp.info.itags(1:2);
    for i = 1 : n - 1
        cnPsi(i) = contract(signOp, 2, psi(i), 2, [2 1 3]);
        cnPsi(i).info.itags{2} = strcat(int2str(i), 's');
    end
    fn = sqrt(2) * S(2);
    if (n ~= length(psi))
        fn.info.itags = {strcat(int2str(n), 's'), strcat(int2str(n), 's*'), strcat(int2str(n), 'm*')};
        if (~isempty(cnPsi(n)))
            cnPsi(n).info.itags{3} = strcat(int2str(n), 'a', cnPsi(n).info.itags{3});
            cnPsi(n) = contract( ...
                contract(fn, 2, psi(n), 2, [3 1 2 4]), '34', ...
                getIdentity(fn, 3, psi(n), 3), '12');
            cnPsi(n).info.itags{3} = strcat(int2str(n), 'a', cnPsi(n).info.itags{3});
        end
    else
        fn.Q{3} = -1 * fn.Q{3};
        fn.info.itags = {strcat(int2str(n), 's'), strcat(int2str(n), 's*'), strcat(int2str(n), 'm')};
        cnPsi(n) = contract( ...
            contract(fn, 2, psi(n), 2, [3 1 2 4]), '34', ...
            getIdentity(fn, 3, psi(n), 3), '12*');
        if (~isempty(cnPsi(n)))
            cnPsi(n).info.itags{3} = strcat(int2str(n), 'a', cnPsi(n).info.itags{3});
        end
    end
   % We changed Q by -2 when applying S(2) on psi(n), we now need to update Qs in the rest of psi.
    for i = n+1 : length(psi)
        cnPsi(i) = psi(i);
        cnPsi(i).Q{1} = cnPsi(i).Q{1} - 2;
        cnPsi(i).Q{3} = cnPsi(i).Q{3} - 2;
    end 
end