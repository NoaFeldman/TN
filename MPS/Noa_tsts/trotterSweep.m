function psi = trotterSweep(psi, dtReal, dtIm, H, orthoOpts)
    % 
    trotterGates = getTrotterGates(H, dtReal, dtIm);
    for k = length(psi) - 1 : -1 : 1
        psi = applyHPair(trotterGates, k, psi, '<<', orthoOpts);
    end
    for k = 1 : length(psi) - 1
        psi = applyHPair(trotterGates, k, psi, '>>', orthoOpts);
    end
end
    