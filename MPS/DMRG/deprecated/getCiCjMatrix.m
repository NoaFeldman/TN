function M = getCiCjMatrix(psi, L)
    % Creates M_ij = <psi|Ci^dagger Cj|psi>
    % M will be an L*L matrix of the first L sites of the chain.
    M = zeros(L, L);
    workPsi = psi;
    canonCenter = length(psi);
    for i = 1 : L
        while canonCenter > i
            [workPsi, canonCenter] = moveCanonicationCenter(workPsi, canonCenter, '<<');
        end
        CiState = getCnForSiteI(workPsi, i, i);
        M(i, i) = getscalar(contract(CiState, '1234', CiState, '1234*'));
        if (i < L)
            [workPsi, canonCenter] = moveCanonicationCenter(workPsi, canonCenter, '>>');
        end
        % TODO ugly, fix it! I only do it because I want to left canonicize
        % the open site, so we have the 'm' leg and the site needs to be
        % contracted by different indices.
        CiState = getCnForSiteI(workPsi, i, i);
        % open = 
        %   ___
        %  | i |--a*
        %    |
        %  | Ci|--m
        %    |
        %  | Cj|
        %    |
        %  |_i_|--a
        % Which will grow by one column n each step.
        open = contract(CiState, '12', getCnForSiteI(workPsi, i+1, i), '12*');
        for j = i+1 : L
            CjState = getCnForSiteI(workPsi, j, j);
            temp = contract(contract(workPsi(j), 3, getIdentity(workPsi(j), 3), '1*'), '23', ...
                CjState, '23*');
            M(j, i) = getscalar(contract(open, '123', temp, '132'));
            M(i, j) = conj(M(j, i));
            
            if (j < L)
                [workPsi, canonCenter] = moveCanonicationCenter(workPsi, canonCenter, '>>');
            end
            open = contract(contract(open, 1, workPsi(j), 1), '23', ...
                            getCnForSiteI(workPsi, j+1, j), '12*', [2 1 3]);
        end
     end
end

function [psi, center] = moveCanonicationCenter(psi, l, dir)
% If dir == '<<', do:
%      _________          _______
%  ->-|psi(l -1)|->-  ->-|psi(l)|-<-   Goes to
%          |                 |
%      ________             
%  ->-|____M___|-<-   Goes to 
%       |    |
%
%      ________          ______
%  ->-|psi(l-1)|-<-  -<-|psi(l)|-<-
%          |                |
    if (dir == '<<')
        opts = {'stol', 1e-7};
        M = contract(psi(l-1), 3, psi(l), 1);
        psi = decomposeAndTruncate(M, l-1, psi, dir, opts);
        center = l-1;
    end
    % equivalent for dir = '>>'
    if (dir == '>>')
        opts = {'stol', 1e-7};
        M = contract(psi(l), 3, psi(l+1), 1);
        psi = decomposeAndTruncate(M, l, psi, dir, opts);
        center = l+1;    
    end
end

function cnPsi = getCnForSiteI(psi, n, i)
    [S,IS]=getLocalSpace('Spin',0.5,'-A');
    if (i < n)
        signOp = -2 * S(1);
        signOp.Q = signOp.Q(1:2);
        signOp.info.itags = signOp.info.itags(1:2);
        cnPsi = contract(signOp, 2, psi(i), 2, [2 1 3]);
        cnPsi.info.itags{2} = strcat(int2str(i), 's');
    elseif (i == n)
        s = sqrt(2) * S(2);
        s.info.itags = {strcat(int2str(n), 's'), strcat(int2str(n), 's*'), 'm*'};
        cnPsi =  contract(s, 2, psi(n), 2, [3 1 4 2]);
    else
        cnPsi = psi(i);
    end    
end