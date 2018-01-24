function [psi, H, HR, HL] = myStartup(N, h, JPM, JZ)
    path(path, [pwd, '/MPSPACK_v3.0']);
    path(path, [pwd, '/Noa_tsts']);
    startup;
    psi = getStartupState(N);
%     % Get the N-1 site to be the working site
%     [psi, workingSiteIndex] = shiftWorkingSite(psi, N, '<<');
    H = getH(N, h, JPM, JZ);
    HR(N) = getHLR(H, psi, N+1, '<<');
%    HR(N-1)   = getHLR(H, psi, N, '<<', HR(N));
    HL(1) = getHLR(H, psi, 0, '>>');
    for l = 1 : N-2
        HL(l+1) = getHLR(H, psi, l, '>>', HL(l));
    end