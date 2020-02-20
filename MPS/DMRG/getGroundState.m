function [psi, H, HR, HL, EGS, Es] = getGroundState(N, params, d, m, bc, psiInitial)
    % Allow for a different startup state (for example for finding first
    % excited state)
    if nargin == 9
        [psi, H, HR, HL, HR2, HL2] = myStartup(N, params, d, psiInitial, m, bc);
    else
        [psi, H, HR, HL, HR2, HL2] = myStartup(N, params, d);
    end
    % Find ground state
    ECurr = 0;
    EError = 1e-6; %1e-8;
    maxBondDimension = 256;
    opts = {'Nkeep', 32};
    for i=1:500
        EForm = ECurr;
        [HL, HR, HL2, HR2, psi, ~] = dmrgSweep(HL, HR, HL2, HR2, H, psi, '<<', opts);
        [HL, HR, HL2, HR2, psi, ECurr] = dmrgSweep(HL, HR, HL2, HR2, H, psi, '>>', opts);
        Es(i) = ECurr;
        disp(i)
        if stepConverged(ECurr, EForm, EError)
            break;
        end
        if mod(i, 10) == 0
            opts{2} = min(opts{2}*2, maxBondDimension);
        end
        if (i == 500)
            disp(['Sweeped 500 times and still not converged, ECurr = ' num2str(ECurr) ...
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