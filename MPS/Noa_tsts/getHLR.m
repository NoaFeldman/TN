function HLR = getHLR(H, psi, l, dir, HLR)
    % Returns <H> for the lth site:
    % If l is in the begining of the chain, returns
    %  _
    % | |--
    % | |
    % | |-- for dir='>>', and the miror QSpace for dir = '<<'
    % | |
    % |_|--
    %  
    % else, performs
    %  _        _
    % | |--  --| |--
    % | |      | |  
    % | |--  --| |--
    % | |      | |
    % |_|--  --|_|--
    %
    % HLR.identityChain is I x I x I ... for all sites contracted (two degree tensor)
    % HLR.opSum  is H(1).single x I x I ... + I x H(2).single x I ... +
    %   H.l2r(1) x H.r2l(2) x I ...  (two degree tensor)
    % HLR.openOp is I x I x ... x H(l).l2r  (three degree tensor)
    if (strcmp(dir, '>>'))
        if (l == 0)
            HLR.identityChain = getIdentity(psi(1), 1, [2, 1]);
            HLR.opSum = QSpace;
            HLR.openOp = QSpace;
            return;
        end
        % TODO document with drawings
        thisSiteHSingle = contract(contract(H.single(l), 1, psi(l), 2), 1, psi(l), '2*');
        thisSiteHSingle = contract(HLR.identityChain, '12', thisSiteHSingle, '13');
        thisSiteIdentity = contract(contract(H.identity(l), 1, psi(l), 2), 1, psi(l), '2*');
        if (l ~= 1)
            HLR.opSum = contract(HLR.opSum, '12', thisSiteIdentity, '13');
            thisSiteHR2L = contract(contract(H.r2l(l), 1, psi(l), 2), 1, psi(l), '2*');
            HLR.opSum  = HLR.opSum + ...
                contract(HLR.openOp, '123', thisSiteHR2L, '124');
        end
        thisSiteHL2R = contract(contract(H.l2r(l), 1, psi(l), 2), 1, psi(l), '2*');
        HLR.openOp = contract(HLR.identityChain , '12', thisSiteHL2R, '24');
        
        HLR.opSum = HLR.opSum + thisSiteHSingle;
        HLR.identityChain = contract(HLR.identityChain, '12', thisSiteIdentity, '13'); 
    end
    if (strcmp(dir, '<<'))
        if (l == length(psi) + 1)
            HLR.identityChain = getIdentity(psi(length(psi)), 3);
            HLR.opSum = QSpace;
            HLR.openOp = QSpace;
            return;
        end
        % TODO document with drawings
        thisSiteHSingle = contract(contract(H.single(l), 1, psi(l), 2), 1, psi(l), '2*');
        thisSiteHSingle = contract(HLR.identityChain, '12', thisSiteHSingle, '42');
        thisSiteIdentity = contract(contract(H.identity(l), 1, psi(l), 2), 1, psi(l), '2*');
        if (l ~= length(psi))
            HLR.opSum = contract(HLR.opSum, '12', thisSiteIdentity, '42');
            thisSiteHL2R = contract(contract(H.l2r(l), 1, psi(l), 2), 1, psi(l), '2*');
            HLR.opSum  = HLR.opSum + ...
                contract(HLR.openOp, '123', thisSiteHL2R, '135');
        end
        thisSiteHR2L = contract(contract(H.r2l(l), 1, psi(l), 2), 1, psi(l), '2*');
        HLR.openOp = contract(HLR.identityChain , '12', thisSiteHR2L, '53');
        HLR.opSum = HLR.opSum + thisSiteHSingle;
        HLR.identityChain = contract(HLR.identityChain, '12', thisSiteIdentity, '42'); 
    end
    
    