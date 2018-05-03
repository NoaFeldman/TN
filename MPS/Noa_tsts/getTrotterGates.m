function trotterGates = getTrotterGates(H, dtReal, dtIm)
    % returns the building blocks for a trotter step.
    % pairOp: 
    %     s_k,k+1
    %   _____|_____
    %  |    H_k    |
    %  |___________|
    %        |
    %    s_k,k+1
    % 
    % IDPair:
    %   s_k     s_k+1
    %   _|_______|_
    %  |  I_k,k+1  |
    %  |___________|
    %        |
    %      s_k,k+1
    % 
    N = length(H.single);
    trotterGates = QSpace(N - 1);
    for k = 1 : N - 1
        IDPair = getIdentity(H.identity(k), 1, H.identity(k+1), 1);
        IDPair.info.itags{3} = ...
            strcat(int2str(k), int2str(k+1), 's', IDPair.info.itags{3});
        pairOp = contract(H.l2r(k), 3, H.r2l(k+1), 3);
        trotterGates(k) = contract(pairOp, '24', IDPair, '12');
        trotterGates(k) = trotterGates(k) + ...
            contract(H.single(k), 2, IDPair, 1);
        if (k == N - 1)
        trotterGates(k) = trotterGates(k) + ...
            contract(H.single(k+1), 2, IDPair, 2, [2 1 3]);
        end
        trotterGates(k) = contract(trotterGates(k), '12', IDPair, '12*');
        IDPairKronned = contract(IDPair, '12', IDPair, '12*');
        trotterGates(k) = trotterGates(k) + IDPairKronned * 1e-25;
         % exponentiate
        for l=1:length(trotterGates(k).data)
            trotterGates(k).data{l} = ...
                expm(-0.5j * complex(dtReal, dtIm).* trotterGates(k).data{l});            
        end
        % bring to 4 rank tensor form again
        trotterGates(k) = contract( ...
            contract(IDPair, 3, trotterGates(k), 2), 3, IDPair, '3*', [3 4 1 2]);
    end
end