function [HL, HR, psi] = dmrgSweep(HL, HR, H, psi, dir, Nkeep)
    if (strcmp(dir, '>>'))
        k = 1;
        while(k < length(psi))
            [HL, HR, psi, k] = dmrgStep(HL, HR, H, psi, k, dir, Nkeep);
        end
    else
        k = length(psi);
        while(k > 1)
            [HL, HR, psi, k] = dmrgStep(HL, HR, H, psi, k, dir, Nkeep);
        end
    end

            