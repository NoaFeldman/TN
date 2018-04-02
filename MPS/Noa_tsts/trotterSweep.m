function [psi, truncErr] = trotterSweep(psi, dtReal, dtIm, H, opts)
    % 
    truncErr = 0;
%     trotterGates = getTrotterGates(H, dtReal, dtIm);
%     for k = length(psi) - 1 : -1 : 1
%         [err, psi] = applyHPair(trotterGates, k, psi, '<<', opts);
%         if (err > truncErr)
%             truncErr = err;
%         end
%     end
%     for k = 1 : length(psi) - 1
%         [err, psi] = applyHPair(trotterGates, k, psi, '>>', opts);
%         if (err > truncErr)
%             truncErr = err;
%         end
%     end
%     TODO just a test to see what causes <N_A> shift from 0. 
    workingSiteIndex = length(psi);
    while workingSiteIndex > 1
        [psi, workingSiteIndex] = shiftWorkingSite(psi, workingSiteIndex, '<<');
    end
    trotterGates = getTrotterGates(H, dtReal, dtIm);
    for k = 1 : length(psi) - 1
        [err, psi] = applyHPair(trotterGates, k, psi, '>>', opts);
        if (err > truncErr)
            truncErr = err;
        end
    end
    for k = length(psi) - 1 : -1 : 1
        [err, psi] = applyHPair(trotterGates, k, psi, '<<', opts);
        if (err > truncErr)
            truncErr = err;
        end
    end
    while workingSiteIndex < length(psi) 
        [psi, workingSiteIndex] = shiftWorkingSite(psi, workingSiteIndex, '>>');
    end
end