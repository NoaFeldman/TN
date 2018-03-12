function [psi, H, HR, HL] = myStartup(N, h, JPM, JZ, m)
    path(path, [pwd, '/MPSPACK_v3.0']);
    path(path, [pwd, '/Noa_tsts']);
    startup;
    psi = getStartupState(N, m);
    H = getH(N, h, JPM, JZ);
    HR(N) = getHLR(H, psi, N+1, '<<');
    HL(1) = getHLR(H, psi, 0, '>>');
    for l = 1 : N-2
        HL(l+1) = getHLR(H, psi, l, '>>', HL(l));
    end
end