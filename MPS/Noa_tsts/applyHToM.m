function Hv = applyHToM(HR, HL, H, M, k)
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
    Hv = Hv + contract(M, 4, HR.opSum, 1);
    Hv = Hv + contract(HL.opSum, 1, M, 1);

    % Add I(Left) x h.single(k1) x h.identity(k2) x
    %   I(Right)
    % And I(Left) x h.identity(k1) x h.single(k2) x
    %   I(Right)
    Hv = Hv + contract(M, 2, H.single(k1), 2, [1 4 2 3]);
    Hv = Hv + contract(M, 3, H.single(k2), 2, [1 2 4 3]);
    % Add HL.openOp x h.r2l(k1) x h.identity(k2) x I(Right)
    % And I(Left) x h.identity(k1) x h.l2r(k2) x HR.openOp
    HK1R2L = contract(H.r2l(k1), 2, M, 2, [2 3 1 4 5]);
    Hv = Hv + contract(HL.openOp, '12', HK1R2L, '12');
    HK2L2R = contract(H.l2r(k2), 2, M, 3, [3 4 1 2 5]);
    Hv = Hv + contract(HK2L2R, '45', HR.openOp, '12');
    % Add I(Left) x h.l2r(k1) x h.r2l(k2) x I(Right)
    HK1K2 = contract(M, 2, H.l2r(k1), 2, [1 4 5 2 3]);
    Hv = Hv + contract(HK1K2, '34', H.r2l(k2), '32', [1 2 4 3]);