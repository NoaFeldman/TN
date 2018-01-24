function exp = getExpectationValue(psi, op, k)
    % Calculate <op>.
    % k is the site which is mixed canonical (must be 1 or N).
    % We construct < psi | op psi > similarly to the way we construct
    % HLR.
    
    if (k == length(psi))
        opSumL = QSpace;
        for i = 1 : k
            opPsi = contract(psi(i), 2, op(i), 2);
            opSumL = contract(contract(opSumL, 1, psi(i), 1), '12', psi(i), '12*');
            opSumL = opSumL + contract(opPsi, '13', psi(i), '12*');
        end
        exp = getscalar(contract(opSumL, '12', getIdentity(opSumL, 1), '21'));
    else
        % TODO
    end
        
        