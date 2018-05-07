function [psi, workingSiteIndex] = shiftWorkingSite(psi, workingSiteIndex, dir)
    % canonicize psi(workingSiteIndex) according to dir, contract the
    % remaining to the next site and promote workingSiteIndex by (1 * dir).
    % Assuming there is a site at (workingSiteIndex + 1 * dir).
    
    opts = {'Nkeep', 1024, 'stol', 1e-16};
    if (strcmp(dir, '>>'))
        [r, psi(workingSiteIndex)] = orthoQS(psi(workingSiteIndex), [1, 2], '>>', opts{:});
        workingSiteIndex = workingSiteIndex + 1;
        psi(workingSiteIndex) = contract(QSpace(r), '1', psi(workingSiteIndex), '1');        
    else if (strcmp(dir, '<<'))
        [psi(workingSiteIndex), l] = orthoQS(psi(workingSiteIndex), 1, '<<', opts{:});
        workingSiteIndex = workingSiteIndex - 1;
        psi(workingSiteIndex) = contract(psi(workingSiteIndex), '3', QSpace(l), '2');
        end
    end
    