function [psi, H, HR, HL] = getGroundState(N, h, JPM, JZ, m, opts)
    [psi, H, HR, HL] = myStartup(N, h, JPM, JZ, m);
    % Find ground state
    ECurr = 0;
    for i=1:100
        disp(strcat('sweep num = ', int2str(i)));
        psi(length(psi) - 1)
        psi(length(psi))
        EForm = ECurr;
        [HL, HR, psi, ~] = dmrgSweep(HL, HR, H, psi, '<<', opts);
        [HL, HR, psi, ECurr] = dmrgSweep(HL, HR, H, psi, '>>', opts);
        if (abs(ECurr - EForm)/abs(ECurr) < 1e-5)
            break;
        end
        if (i == 100)
            disp(['Sweeped 100 times and still not converged, Delta = ' ...
                num2str(Delta) ', ECurr = ' num2str(ECurr) ...
                ', EForm = ' num2str(EForm)]);
        end
    end
end