function tDependentNegs(inDir, inInds, outDir)
    path(path, ['../MPSPACK_v3.0']);
    mkdir(outDir);
    for i = inInds
        d = load(strcat(inDir, '/psiAtStep', int2str(i)));
        currPsi = shrink(d.psi);
        l = 16;
        saveNegativitySpectrum(currPsi, l, 'symm', strcat(outDir, '/step', int2str(i)));
    end
end

function shrinked = shrink(currPsi)
    k = length(currPsi);
    while k > 1
        [currPsi, k, ~] = shiftWorkingSite(currPsi, k, '<<', {'Nkeep', 20});
    end
    while k < length(currPsi)
        [currPsi, k, ~] = shiftWorkingSite(currPsi, k, '>>', {'Nkeep', 20});
    end
    shrinked = currPsi;
end
    
    