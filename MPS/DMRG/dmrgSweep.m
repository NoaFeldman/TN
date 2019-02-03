function [HL, HR, HL2, HR2, psi, E] = dmrgSweep(HL, HR, HL2, HR2, H, psi, dir, opts)
    if (strcmp(dir, '>>'))
        k = 1;
        while(k < length(psi))
            [HL, HR, HL2, HR2, psi, E, k] = dmrgStep(HL, HR, HL2, HR2, H, psi, k, dir, opts);
        end
    else
        k = length(psi);
        while(k > 1)
            [HL, HR, HL2, HR2, psi, E, k] = dmrgStep(HL, HR, HL2, HR2, H, psi, k, dir, opts);
        end
    end
end

    
    
    
    
    
    
    
    
          