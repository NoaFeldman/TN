function Hv = applyH(HR, HL, H, psi, k, dir)
    % HL, HR represent the two blocks of <H> from both sides, H(k) and
    % H(k+dir) the MPO of the kth and (k+dir)th sites, and psi(k) + psi(k+dir)
    % are the matrices to be optimized.
    % If dir == '>>' then the second MPS is on the left (i.e. k, k+1),
    % otherwise we use k, k-1. 
    %  __                       __
    % |  |---               ---|  |
    % |HL|    _|__   __|___    |HR|
    % |  |---|H(k)|-|H(k+1)|---|  |
    % |  |     |       |       |  |
    % |  |                     |  |
    % |  |  ___|___  __|_____  |  |
    % |__|-|psi(k)|-|psi(k+1)|-|__|
    % The tensor returned is of degree 4 as seen in the drawing.
    Hv = QSpace;
    k1 = k;
    if (strcmp(dir, '>>'))
        k2 = k+1;
    else
        k2 = k - 1;
    end
    
    k1SiteIdentity = contract(H.identity(k1), 2, psi(k1), '2*');
    k2SiteIdentity = contract(H.identity(k2), 2, psi(k2), '2*');    
    % Add HL.opSum x h.identity(k1) x h.identity(k2) x
    %   HR.identityChain
    % and HL.identityChain x h.identity(k1) x h.identity(k2) x
    %   HR.opSum
    Hv = Hv + ...
        contract(contract(contract(HL.identityChain, 2, k1SiteIdentity, 2), 3, ...
            k2SiteIdentity, 2), 4, HR.opSum, 2);
    Hv = Hv + ...
        contract(contract(contract(HL.opSum, 2, k1SiteIdentity, 2), 3, ...
            k2SiteIdentity, 2), 4, HR.identityChain, 2);
    % Add HL.identityChain x h.single(k1) x h.identity(k2) x
    %   HR.identityChain
    % And HL.identityChain x h.identity(k1) x h.single(k2) x
    %   HR.identityChain
    HK1SiteSingle  = contract(H.single(k1), 2, psi(k1), '2*');    
    HK1SiteSingle  = contract(HK1SiteSingle, 3, k2SiteIdentity, 2);
    HK1SiteSingle  = contract(HL.identityChain, 2, contract(HK1SiteSingle, 4, HR.identityChain, 2), 2);
    Hv = Hv + HK1SiteSingle;
    HK2SiteSingle  = contract(H.single(k2), 2, psi(k2), '2*');
    HK2SiteSingle  = contract(k1SiteIdentity, 3, HK2SiteSingle, 2);
    HK2SiteSingle  = contract(HL.identityChain, 2, contract(HK2SiteSingle, 4, HR.identityChain, 2), 2);
    Hv = Hv + HK2SiteSingle;
    % TODO suitable only for dir = '>>'
    % Add HL.openOp x h.r2l(k1) x h.identity(k2) x HR.identityChain
    % And HL.identityChain x h.identity(k1) x h.l2r(k2) x HR.openOp
    HK1R2L = contract(H.r2l(k1), 2, psi(k1), '2*');
    HK1R2L = contract(HL.openOp, '13', HK1R2L, '23');
    HK1R2L = contract(HK1R2L, 3, k2SiteIdentity, 2);
    HK1R2L = contract(HK1R2L, 4, HR.identityChain, 2);
    Hv = Hv + HK1R2L;
    HK2L2R = contract(H.l2r(k2), 2, psi(k2), '2*');
    HK2L2R = contract(HK2L2R, '24', HR.openOp, '13');
    HK2L2R = contract(k1SiteIdentity, 3, HK2L2R, 2);
    HK2L2R = contract(HL.identityChain, 2, HK2L2R, 2);
    Hv = Hv + HK2L2R;
    % Add HL.identityChain x h.l2r(k1) x h.r2l(k2) x HR.identityChain
    HK1L2R = contract(H.l2r(k1), 2, psi(k1), '2*');
    HK2R2L = contract(H.r2l(k2), 2, psi(k2), '2*');
    HK1K2 = contract(HK1L2R, '24', HK2R2L, '23');
    HK1K2 = contract(contract(HL.identityChain, 2, HK1K2, 2), 4, HR.identityChain, 2);
    Hv = Hv + HK1K2;