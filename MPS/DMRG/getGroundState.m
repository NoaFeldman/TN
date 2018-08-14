function [psi, H, HR, HL] = getGroundState(N, h, JPM, JZ, m)
    [psi, H, HR, HL] = myStartup(N, h, JPM, JZ, m);
    % Find ground state
    ECurr = 0;
    EError = 1e-6;
    maxBondDimension = 128;
    opts = {'Nkeep', 4};
    for i=1:400
        EForm = ECurr;
        [HL, HR, psi, ~] = dmrgSweep(HL, HR, H, psi, '<<', opts);
        [HL, HR, psi, ECurr] = dmrgSweep(HL, HR, H, psi, '>>', opts);
        if (stepConverged(ECurr, EForm, EError))
            [HL, HR, psi, dmrgConverged] = isConverged(HL, HR, H, psi, opts, ECurr, EError);
            if (dmrgConverged)
                break;
            end
            opts{2} = opts{2}*2;
        end        
        if (i == 400)
            disp(['Sweeped 100 times and still not converged, ECurr = ' num2str(ECurr) ...
                ', EForm = ' num2str(EForm) ...
                ', NKeep = ' num2str(opts{2})]);
        end
    end
end

function converged = stepConverged(ECurr, EForm, EError)
    converged = abs(ECurr - EForm)/abs(ECurr) < EError;
end

function [HL, HR, psi, dmrgConverged] = isConverged(HL, HR, H, psi, opts, ECurr, EError)
    EForm = ECurr;
    opts{2} = 2 * opts{2};
    [HL, HR, psi, ~] = dmrgSweep(HL, HR, H, psi, '<<', opts);
    [HL, HR, psi, ECurr] = dmrgSweep(HL, HR, H, psi, '>>', opts);
    if (stepConverged(ECurr, EForm, EError))
        dmrgConverged = true;
    else
        dmrgConverged = false;
    end      
end