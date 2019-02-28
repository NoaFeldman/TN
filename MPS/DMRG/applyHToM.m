function Hv = applyHToM(HR, HL, HR2, HL2, H, M, k, bc)
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
    %   I(Right)
    % and I(Left) x h.identity(k1) x h.identity(k2) x
    %   HR.opSum
    if strcmp(bc, 'open')
        Hv = Hv + contract(M, 4, HR.opSum, 1);
        Hv = Hv + contract(HL.opSum, 1, M, 1);
    else
        Hv = Hv + contract(M, 6, HR.opSum, 1);
        Hv = Hv + contract(HL.opSum, 1, M, 1);        
    end

    % Add I(Left) x h.single(k1) x h.identity(k2) x
    %   I(Right)
    % And I(Left) x h.identity(k1) x h.single(k2) x
    %   I(Right)
    if strcmp(bc, 'open')
        Hv = Hv + contract(M, 2, H.single(k1), 2, [1 4 2 3]);
        Hv = Hv + contract(M, 3, H.single(k2), 2, [1 2 4 3]);
    else
        Hv = Hv + contract(M, '23', H.single(k1), '24', [1 5 6 2 3 4]);
        Hv = Hv + contract(M, '45', H.single(k2), '24', [1 2 3 5 6 4]);        
    end
    % Add HL.openOp x h.r2l(k1) x h.identity(k2) x I(Right)
    % And I(Left) x h.identity(k1) x h.l2r(k2) x HR.openOp
    if strcmp(bc, 'open')
        HK1R2L = contract(H.r2l(k1), 2, M, 2, [2 3 1 4 5]);
        Hv = Hv + contract(HL.openOp, '12', HK1R2L, '12');
        HK1R2L2 = contract(H.r2l2(k1), 2, M, 2, [2 3 1 4 5]);
        Hv = Hv + contract(HL2.toClose, '12', HK1R2L2, '12');
        HK2L2R = contract(H.l2r(k2), 2, M, 3, [3 4 1 2 5]);
        Hv = Hv + contract(HK2L2R, '45', HR.openOp, '12');
        HK2L2R2 = contract(H.l2r2(k2), 2, M, 3, [3 4 1 2 5]);
        Hv = Hv + contract(HK2L2R2, '45', HR2.toClose, '12');
    else
        HK1R2L = contract(H.r2lUp(k1), '24', M, '23', [3 4 1 2 5 6 7]);
        Hv = Hv + contract(HL.openOpUp, '12', HK1R2L, '12');
        HK1R2L = contract(H.r2lDown(k1), '24', M, '23', [3 4 1 2 5 6 7]);
        Hv = Hv + contract(HL.openOpDown, '12', HK1R2L, '12');
        
        HK1R2L2 = contract(H.r2l2Up(k1), '24', M, '23', [3 4 1 2 5 6 7]);
        Hv = Hv + contract(HL2.toCloseUp, '12', HK1R2L2, '12');
        HK1R2L2 = contract(H.r2l2Down(k1), '24', M, '23', [3 4 1 2 5 6 7]);
        Hv = Hv + contract(HL2.toCloseDown, '12', HK1R2L2, '12');
        
        HK2L2R = contract(H.l2rUp(k2), '24', M, '45', [4 5 6 1 2 3 7]);
        Hv = Hv + contract(HK2L2R, '67', HR.openOpUp, '12');        
        HK2L2R = contract(H.l2rDown(k2), '24', M, '45', [4 5 6 1 2 3 7]);
        Hv = Hv + contract(HK2L2R, '67', HR.openOpDown, '12');        
        
        HK2L2R2 = contract(H.l2r2Up(k2), '24', M, '45', [4 5 6 1 2 3 7]);
        Hv = Hv + contract(HK2L2R2, '67', HR2.toCloseUp, '12');        
        HK2L2R2 = contract(H.l2r2Down(k2), '24', M, '45', [4 5 6 1 2 3 7]);
        Hv = Hv + contract(HK2L2R2, '67', HR2.toCloseDown, '12');        
    end
    % Add I(Left) x h.l2r(k1) x h.r2l(k2) x I(Right)
    if strcmp(bc, 'open')
        HK1K2 = contract(M, 2, H.l2r(k1), 2, [1 4 5 2 3]);
        Hv = Hv + contract(HK1K2, '34', H.r2l(k2), '32', [1 2 4 3]);
    else
        HK1K2 = contract(M, '23', H.l2rUp(k1), '24', [1 5 6 2 3 4 7]);
        Hv = Hv + contract(HK1K2, '457', H.r2lUp(k2), '245', [1 2 3 5 6 4]);
        HK1K2 = contract(M, '23', H.l2rDown(k1), '24', [1 5 6 2 3 4 7]);
        Hv = Hv + contract(HK1K2, '457', H.r2lDown(k2), '245', [1 2 3 5 6 4]);
    end
    % Add cotribution of HL2.toContinue with right site of M and its mirror
    if strcmp(bc, 'open')
        temp = contract(HL2.toContinue, 2, M, 1);
        Hv = Hv  + contract(temp, '14', H.r2l2(k2), '32', [1 2 4 3]);
        temp = contract(HR2.toContinue, 2, M, 4);
        Hv = Hv  + contract(temp, '14', H.l2r2(k1), '32', [2 4 3 1]);
    else
        temp = contract(HL2.toContinueUp, 2, M, 1);
        Hv = Hv  + contract(temp, '156', H.r2l2Up(k2), '524', [1 2 3 4 6 7 5]);
        temp = contract(HL2.toContinueDown, 2, M, 1);
        Hv = Hv  + contract(temp, '156', H.r2l2Down(k2), '524', [1 2 3 4 6 7 5]);
        temp = contract(HR2.toContinueUp, 2, M, 6);
        Hv = Hv  + contract(temp, '145', H.l2r2Up(k1), '524', [2 5 6 3 4 1]);
        temp = contract(HR2.toContinueDown, 2, M, 6);
        Hv = Hv  + contract(temp, '145', H.l2r2Down(k1), '524', [2 5 6 3 4 1]);
    end
end