function [psi, H, HR, HL, HR2, HL2] = myStartup(N, params, d, psi, m, bc)
    path(path, [pwd, '/MPSPACK_v3.0']);
    path(path, [pwd, '/DMRG']);
    startup;
    if nargin == 3
        m = 0;
        bc = 'open';
        psi = getStartupState(N, d, m, bc);
    end
    if strcmp(bc, 'periodic')
        N = N/2;
    end
    H = getH(N, d, params, bc);
    [HR(N), HR2(N)] = getHLR(H, psi, N+1, '<<');
    [HL(1), HL2(1)] = getHLR(H, psi, 0, '>>');
    for l = 1 : N-2
        [HL(l+1), HL2(l+1)] = getHLR(H, psi, l, '>>', HL(l), HL2(l));
    end
end