function [psi, H, HR, HL] = getGroundState(N, h, JPM, JZ, m, opts)
    [psi, H, HR, HL] = myStartup(N, h, JPM, JZ, m);
    % Find ground state
    ECurr = 0;
    for i=1:100
        EForm = ECurr;
        [HL, HR, psi, ~] = dmrgSweep(HL, HR, H, psi, '<<', opts);
        [HL, HR, psi, ECurr] = dmrgSweep(HL, HR, H, psi, '>>', opts);
        if (abs(ECurr - EForm)/abs(ECurr) < 1e-7)
            break;
        end        
        if (i == 100)
            disp(['Sweeped 100 times and still not converged, ECurr = ' num2str(ECurr) ...
                ', EForm = ' num2str(EForm)]);
        end
    end
end