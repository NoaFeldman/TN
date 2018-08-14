function groundStateSpectrumDMRG(Ls, JZs)
    for L = Ls
        for jz = JZs
            [psi, ~, ~, ~] = getGroundState(L, 0, 1, jz, 0);
            saveRDMSpectrum(strcat('groundStateSpectrumDMRG_', int2str(L), '_', int2str(L/2), '_0', int2str(10*jz)), psi, L/2);
        end
    end
end