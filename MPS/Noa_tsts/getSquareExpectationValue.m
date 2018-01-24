function exp = getSquareExpectationValue(psi, op, k)
    % Calculate expectation value of op^2.
    % k is the site which is mixed canonical (must be 1 or N).
    % We construct <op psi | op psi > similarly to the way we construct
    % HLR.
    
    % We keep 8 open tensors from each side:
    % upSum - we have operated a term in op on some site on the upper state
    % (psi), but not on the lower site (psi^dagger).
    % downSum - analogous
    % closedSum - we have operated on both sites
    % openOpUp - 
    % openOpDown - 
    % openUpclosedDown -
    % openDownClosedUp -
    % openBoth - 
    if (k == length(psi))
        upSumL = QSpace;
        downSumL = QSpace;
        closedSumL = QSpace;
        for i = 1 : k
            opPsi = contract(psi(i), 2, op(i), 2);
            closedSumL = contract(contract(closedSumL, 1, psi(i), 1), '12', psi(i), '12*');
            closedSumL = closedSumL + contract(opPsi, '13', opPsi, '13*');
            closedSumL = closedSumL + contract( ...
                contract(upSumL, 1, psi(i), 1), '12', opPsi, '13*');
            closedSumL = closedSumL + contract( ...
                opPsi, '13', contract(downSumL, 2, psi(i), '1*'), '12');

            upSumL = contract( ...
                contract(upSumL, 1, psi(i), 1), '12', psi(i), '12*');
            up = contract(opPsi, '13', psi(i), '12*');
            upSumL = upSumL + up;

            downSumL = contract( ...
                contract(downSumL, 1, psi(i), 1), '12', psi(i), '12*');
            down = contract(psi(i), '12', opPsi, '13*');
            downSumL = downSumL + down;
        end
        exp = getscalar(contract(closedSumL, '12', getIdentity(closedSumL, 1), '12'));
    else
        % TODO
    end
        
        