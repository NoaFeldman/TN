function [psi, truncErr] = trotterSweep(trotterGates, psi, opts)
    % 
    truncErr = 0;
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