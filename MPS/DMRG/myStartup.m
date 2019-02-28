function [psi, H, HR, HL, HR2, HL2] = myStartup(N, h, JPM, JZ, J2PM, J2Z, m, bc)
    path(path, [pwd, '/MPSPACK_v3.0']);
    path(path, [pwd, '/DMRG']);
    startup;
    if strcmp(bc, 'periodic')
        N = N/2;
    end
    psi = getStartupState(N, m, bc);
    H = getH(N, h, JPM, JZ, J2PM, J2Z, bc);
    [HR(N), HR2(N)] = getHLR(H, psi, N+1, '<<');
    [HL(1), HL2(1)] = getHLR(H, psi, 0, '>>');
    for l = 1 : N-2
        [HL(l+1), HL2(l+1)] = getHLR(H, psi, l, '>>', HL(l), HL2(l));
    end
end