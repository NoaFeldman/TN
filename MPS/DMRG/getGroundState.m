function [psi, H, HR, HL, EGS] = getGroundState(N, h, JPM, JZ, J2PM, J2Z, m, bc)
    [psi, H, HR, HL, HR2, HL2] = myStartup(N, h, JPM, JZ, J2PM, J2Z, m, bc);
    % Find ground state
    ECurr = 0;
    EError = 1e-7;
    maxBondDimension = 512;
    opts = {'Nkeep', 4};
    for i=1:1000
        EForm = ECurr;
        [HL, HR, HL2, HR2, psi, ~] = dmrgSweep(HL, HR, HL2, HR2, H, psi, '<<', opts);
        [HL, HR, HL2, HR2, psi, ECurr] = dmrgSweep(HL, HR, HL2, HR2, H, psi, '>>', opts);
        if (stepConverged(ECurr, EForm, EError))
            [HL, HR, HL2, HR2, psi, dmrgConverged] = isConverged(HL, HR, HL2, HR2, H, psi, opts, ECurr, EError);
            if (dmrgConverged)
                break;
            end
            opts{2} = opts{2}*2;
        end        
        if (i == 1000)
            disp(['Sweeped 100 times and still not converged, ECurr = ' num2str(ECurr) ...
                ', EForm = ' num2str(EForm) ...
                ', NKeep = ' num2str(opts{2})]);
        end
    end
    EGS = ECurr;
end

function converged = stepConverged(ECurr, EForm, EError)
    converged = abs(ECurr - EForm)/abs(ECurr) < EError;
end

function [HL, HR, HL2, HR2, psi, dmrgConverged] = isConverged(HL, HR, HL2, HR2, H, psi, opts, ECurr, EError)
    EForm = ECurr;
    opts{2} = 2 * opts{2};
    [HL, HR, HL2, HR2, psi, ~] = dmrgSweep(HL, HR, HL2, HR2, H, psi, '<<', opts);
    [HL, HR, HL2, HR2, psi, ECurr] = dmrgSweep(HL, HR, HL2, HR2, H, psi, '>>', opts);
    if (stepConverged(ECurr, EForm, EError))
        dmrgConverged = true;
    else
        dmrgConverged = false;
    end      
end