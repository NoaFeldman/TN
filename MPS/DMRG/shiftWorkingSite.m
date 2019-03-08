function [psi, workingSiteIndex, truncErr] = shiftWorkingSite(psi, workingSiteIndex, dir, opts)
    % canonicize psi(workingSiteIndex) according to dir, contract the
    % remaining to the next site and promote workingSiteIndex by (1 * dir).
    % Assuming there is a site at (workingSiteIndex + 1 * dir).
    
    if nargin == 3
        opts = {'Nkeep', 1024, 'stol', 1e-16};
    end
    
    if length(psi(1).info.itags) == 3
        bc = 'open';
    else
        bc = 'periodic';
    end

    if (strcmp(dir, '>>'))
        if strcmp(bc, 'open')
            inds = [1 2];
        else
            inds = [1 2 3];
        end
        [r, psi(workingSiteIndex), I] = orthoQS(psi(workingSiteIndex), inds, '>>', opts{:});
        workingSiteIndex = workingSiteIndex + 1;
        psi(workingSiteIndex) = contract(QSpace(r), '1', psi(workingSiteIndex), '1');
        truncErr = I.svd2tr;
    else if (strcmp(dir, '<<'))
        [psi(workingSiteIndex), l, I] = orthoQS(psi(workingSiteIndex), 1, '<<', opts{:});
        workingSiteIndex = workingSiteIndex - 1;
        psi(workingSiteIndex) = contract(psi(workingSiteIndex), '3', QSpace(l), '2');
        truncErr = I.svd2tr;
        end
    end
end
    