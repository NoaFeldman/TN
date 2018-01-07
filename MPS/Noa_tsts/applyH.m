function Hv = applyH(HR, HL, H, M, k)
    % HL, HR represent the two blocks of <H> from both sides, H(k) and
    % H(k+dir) the MPO of the kth and (k+dir)th sites, and psi(k) + psi(k+dir)
    % are the matrices to be optimized.
    % If dir == '>>' then the second MPS is on the left (i.e. k, k+1),
    % otherwise we use k, k-1. 
    %  __   _________________   __
    % |  |-|________M________|-|  |
    % |  |     |       |       |  |
    % |HL|    _|__   __|___    |HR|
    % |  |---|H(k)|-|H(k+1)|---|  |
    % |  |     |       |       |  |
    % |  |                     |  |
    % |__|---               ---|__|
    % The tensor returned is of degree 4 as seen in the drawing.
    Hv = QSpace;
    k1 = k;
    k2 = k+1;
        
    % Add HL.opSum x h.identity(k1) x h.identity(k2) x
    %   HR.identityChain
    % and HL.identityChain x h.identity(k1) x h.identity(k2) x
    %   HR.opSum
    Hv = Hv + ...
        contract(contract(HL.identityChain, 1, M, 1), 4, HR.opSum, 1);
    Hv = Hv + ...
        contract(contract(HL.opSum, 1, M, 1), 4, HR.identityChain, 1);

    % Add HL.identityChain x h.single(k1) x h.identity(k2) x
    %   HR.identityChain
    % And HL.identityChain x h.identity(k1) x h.single(k2) x
    %   HR.identityChain
    HK1SiteSingle  = contract(M, 2, H.single(k1), 1, [1 4 2 3]);    
    Hv = Hv + ...
        contract(contract(HL.identityChain, 1, HK1SiteSingle, 1), ...
            4, HR.identityChain, 1);
    HK2SiteSingle  = contract(M, 3, H.single(k2), 1, [1 2 4 3]);    
    Hv = Hv + ...
        contract(contract(HL.identityChain, 1, HK2SiteSingle, 1), ...
            4, HR.identityChain, 1);
    % Add HL.openOp x h.r2l(k1) x h.identity(k2) x HR.identityChain
    % And HL.identityChain x h.identity(k1) x h.l2r(k2) x HR.openOp
    HK1R2L = contract(H.r2l(k1), 1, M, 2, [2 3 1 4 5]);
    Hv = Hv + ...
        contract(contract(HL.openOp, '12', HK1R2L, '12'), ...
            4, HR.identityChain, 1);
    HK2L2R = contract(H.l2r(k2), 1, M, 3, [3 4 1 2 5]);
    Hv = Hv + ...
        contract(HL.identityChain, 1, ...
        contract(HK2L2R, '45', HR.openOp, '12'), 1);
    % Add HL.identityChain x h.l2r(k1) x h.r2l(k2) x HR.identityChain
    HK1K2 = contract(M, 2, H.l2r(k1), 1, [1 4 5 2 3]);
    HK1K2 =  contract(HK1K2, '34', H.r2l(k2), '31', [1 2 4 3]);
    Hv = Hv + ...
        contract(contract(HL.identityChain, 1, HK1K2, 1), ...
            4, HR.identityChain, 1);