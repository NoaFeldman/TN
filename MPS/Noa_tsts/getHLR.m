function HLR = getHLR(H, psi, l, dir, HLR)
    % TODO
    if (strcmp(dir, '>>'))
        if (l == 1)
            HLR.identities = getIdentity(psi(1), 1, [2, 1]);
            HLR.opSum = QSpace;
            HLR.openOp = QSpace;
        end
        % TODO document with drawings
        thisSiteHSingle = contract(contract(H.single(l), 1, psi(l), 2), 1, psi(l), '2*');
        thisSiteHSingle = contract(HLR.identities, '12', thisSiteHSingle, '13');
        thisSiteIdentity = contract(contract(H.identity(l), 1, psi(l), 2), 1, psi(l), '2*');
        if (l ~= 1)
            HLR.opSum = contract(HLR.opSum, '12', thisSiteIdentity, '13');
            thisSiteHR2L = contract(contract(H.r2l(l), 1, psi(l), 2), 1, psi(l), '2*');
            HLR.opSum  = HLR.opSum + ...
                contract(HLR.openOp, '123', thisSiteHR2L, '124');
        end
        thisSiteHL2R = contract(contract(H.l2r(l), 1, psi(l), 2), 1, psi(l), '2*');
        HLR.openOp = contract(HLR.identities , '12', thisSiteHL2R, '24');
        
        HLR.opSum = HLR.opSum + thisSiteHSingle;
        HLR.identities = contract(HLR.identities, '12', thisSiteIdentity, '13'); 
    end
    if (strcmp(dir, '<<'))
        if (l == length(psi))
            HLR.identities = getIdentity(psi(l), 1, [2, 1]);
            HLR.opSum = QSpace;
            HLR.openOp = QSpace;
        end
        % TODO document with drawings
        thisSiteHSingle = contract(contract(H.single(l), 1, psi(l), 2), 1, psi(l), '2*');
        thisSiteHSingle = contract(HLR.identities, '12', thisSiteHSingle, '13');
        thisSiteIdentity = contract(contract(H.identity(l), 1, psi(l), 2), 1, psi(l), '2*');
        if (l ~= 1)
            HLR.opSum = contract(HLR.opSum, '12', thisSiteIdentity, '13');
            thisSiteHR2L = contract(contract(H.r2l(l), 1, psi(l), 2), 1, psi(l), '2*');
            HLR.opSum  = HLR.opSum + ...
                contract(HLR.openOp, '123', thisSiteHR2L, '124');
        end
        thisSiteHL2R = contract(contract(H.l2r(l), 1, psi(l), 2), 1, psi(l), '2*');
        HLR.openOp = contract(HLR.identities , '12', thisSiteHL2R, '24');
        
        HLR.opSum = HLR.opSum + thisSiteHSingle;
        HLR.identities = contract(HLR.identities, '12', thisSiteIdentity, '13'); 
    end    
    
    