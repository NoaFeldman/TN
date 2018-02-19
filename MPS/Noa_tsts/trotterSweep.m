function [psi, truncErr] = trotterSweep(psi, dtReal, dtIm, H, opts)
    % 
    truncErr = 0;
    trotterGates = getTrotterGates(H, dtReal, dtIm);
    for k = length(psi) - 1 : -1 : 1
        [err, psi] = applyHPair(trotterGates, k, psi, '<<', opts);
        if (err > truncErr)
            truncErr = err;
        end
    end
    for k = 1 : length(psi) - 1
        [err, psi] = applyHPair(trotterGates, k, psi, '>>', opts);
        if (err > truncErr)
            truncErr = err;
        end
    end
end