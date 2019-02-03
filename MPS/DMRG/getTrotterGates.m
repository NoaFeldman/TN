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
    trotterGates.nearest = QSpace(N - 1);
    trotterGates.next = QSpace(N - 2);
    for k = 1 : N - 1
        IDPair = getIdentity(H.identity(k), 1, H.identity(k+1), 1);
        IDPair.info.itags{3} = ...
            strcat(int2str(k), int2str(k+1), 's', IDPair.info.itags{3});
        pairOp = contract(H.l2r(k), 3, H.r2l(k+1), 3);
        trotterGates.nearest(k) = contract(pairOp, '24', IDPair, '12');
        trotterGates.nearest(k) = trotterGates.nearest(k) + ...
            contract(H.single(k), 2, IDPair, 1);
        if (k == N - 1)
        trotterGates.nearest(k) = trotterGates.nearest(k) + ...
            contract(H.single(k+1), 2, IDPair, 2, [2 1 3]);
        end
        trotterGates.nearest(k) = contract(trotterGates.nearest(k), '12', IDPair, '12*');
        IDPairKronned = contract(IDPair, '12', IDPair, '12*');
        trotterGates.nearest(k) = trotterGates.nearest(k) + IDPairKronned * 1e-25;
         % exponentiate
        for l=1:length(trotterGates.nearest(k).data)
            trotterGates.nearest(k).data{l} = ...
                expm(-0.5j * complex(dtReal, dtIm).* trotterGates.nearest(k).data{l});            
        end
        % bring to 4 rank tensor form again
        trotterGates.nearest(k) = contract( ...
            contract(IDPair, 3, trotterGates.nearest(k), 2), 3, IDPair, '3*', [3 4 1 2]);
    end
    
    
    for k = 1 : N - 2
        IDPair = getIdentity(H.identity(k), 1, H.identity(k+2), 1);
        IDPair.info.itags{3} = ...
            strcat(int2str(k), int2str(k+2), 's', IDPair.info.itags{3});
        pairOp = contract(H.l2r2(k), 3, H.r2l2(k+2), 3);
        trotterGates.next(k) = contract(pairOp, '24', IDPair, '12');
        trotterGates.next(k) = trotterGates.next(k) + ...
            contract(H.single(k), 2, IDPair, 1);
        if (k == N - 2)
        trotterGates.next(k) = trotterGates.next(k) + ...
            contract(H.single(k+2), 2, IDPair, 2, [2 1 3]);
        end
        trotterGates.next(k) = contract(trotterGates.next(k), '12', IDPair, '12*');
        IDPairKronned = contract(IDPair, '12', IDPair, '12*');
        trotterGates.next(k) = trotterGates.next(k) + IDPairKronned * 1e-25;
         % exponentiate
        for l=1:length(trotterGates.next(k).data)
            trotterGates.next(k).data{l} = ...
                expm(-0.5j * complex(dtReal, dtIm).* trotterGates.next(k).data{l});            
        end
        % bring to 4 rank tensor form again
        trotterGates.next(k) = contract( ...
            contract(IDPair, 3, trotterGates.next(k), 2), 3, IDPair, '3*', [3 4 1 2]);
    end
end