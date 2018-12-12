function processMPSNeg(inDir, outDir, tSteps, l)
    mkdir outDir;
    for i = 1:length(tSteps)
        t = tSteps(i);
        load(strcat(inDir, '/psiAtStep', int2str(t)));
        saveNegativitySpectrum(psi, l, 'symm', strcat(outDir, '/step', int2str(t)));
    end
end