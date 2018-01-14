function psi = trotterSweep(psi, dtReal, dtIm, H, orthoOpts)
    % 
    HPairs = getHPairs(H, dtReal, dtIm);
    for k = length(psi) - 1 : -1 : 1
        psi = applyHPair(HPairs, k, psi, '<<', orthoOpts);
    end
    for k = 1 : length(psi) - 1
        psi = applyHPair(HPairs, k, psi, '>>', orthoOpts);
    end
end
    