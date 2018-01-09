function [HL, HR, psi, E] = dmrgSweep(HL, HR, H, psi, dir, Nkeep)
    if (strcmp(dir, '>>'))
        k = 1;
        while(k < length(psi))
            [HL, HR, psi, k] = dmrgStep(HL, HR, H, psi, k, dir, Nkeep);
        end
        HRCurr = getHLR(H, psi, length(psi), '^');
        E = getscalar( ...
            contract(getIdentity(HRCurr.opSum, 1), '12', HRCurr.opSum, '12'));
        E = E + getscalar( ...
            contract(getIdentity(HL(length(psi)).opSum, 1, [2 1]), '12', HL(length(psi)).opSum, '12'));
        E = E + getscalar(contract(HRCurr.openOp, '123', HL(length(psi)).openOp, '123'));
    else
        k = length(psi);
        while(k > 1)
            [HL, HR, psi, k] = dmrgStep(HL, HR, H, psi, k, dir, Nkeep);
        end
        HLCurr = getHLR(H, psi, 1, '^');
        E = getscalar( ...
            contract(getIdentity(HLCurr.opSum, 1, [2, 1]), '12', HLCurr.opSum, '12'));
        E = E + getscalar( ...
            contract(getIdentity(HR(1).opSum, 1), '12', HR(1).opSum, '12'));
        E = E + getscalar(contract(HR(1).openOp, '123', HLCurr.openOp, '123'));

    end
end

    
    
    
    
    
    
    
    
          