function processMPSNeg(inDir, outDir, tSteps, l)
    path(path, [pwd '/MPSPACK_v3.0']);
    path(path, [pwd '/DMRG']);
    startup;
    mkdir outDir;
    for i = 1:length(tSteps)
        t = tSteps(i);
        load(strcat(inDir, '/psiAtStep', int2str(t)));
        saveNegativitySpectrum(psi, l, 'symm', strcat(outDir, '/step', int2str(t)));
    end
end