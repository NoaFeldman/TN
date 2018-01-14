function E2 = getHSquared(H, psi, k)
    % Return <H^2>
    % k represents the working site (mixed canonical).
    Hpsi.identityChain = getIdentity(psi(1), 1, [2 1]);
    % opSumUp (opSumDown) is the sum of all combinations with some operator (single or double) 
    % acting on the upper (lower) state, but on the lower (upper) state we have only
    % identities. 
    % opSumBoth is the sum of all combinations with some operator acting on
    % the upper state and some operator acting on the lower.
    Hpsi.opSumUp.curr = QSpace;
    Hpsi.opSumUp.next = QSpace;
    Hpsi.opSumDown.curr = QSpace;
    Hpsi.opSumDown.next = QSpace;
    Hpsi.opSumBoth.curr = QSpace;
    Hpsi.opSumBoth.next = QSpace;
    for i = 1 : k
        identitySite = contract(psi(i), '2', psi(i), '2*');
        pair = contract(psi(i), 3, psi(i+1), 1);
        identityPair = contract(pair, '23', pair, '23*');
        pairOpPsi = ...
            contract(contract(H.l2r(i), 3, H.r2l(i+1), 3), '13', pair, '23', [3 1 2 4]) + ...
        pairOpPsiUp = contract(pairOpPsi, '23', pair, '23*');
        pairOpPsiDown = contract(pair, '23', pairOpPsi, '23*');
        % Add to opSumBoth ID(1 ... i-1) x H(pair)H*(pair) 
        %   + opSumUp x H*(pair) + opSumDown x H(pair)
        Hpsi.opSumBoth.curr = contract(Hpsi.op
        Hpsi.opSumBoth.next = contract(Hpsi.opSumBoth.curr, '12', identityPair, '13');
        Hpsi.opSumBoth.next = Hpsi.opSumBoth + ...
            contract(Hpsi.opSumDown, '12', pairOpPsiUp, '13') + ...
            contract(Hpsi.opSumUp, '12', pairOpPsiDown, '13') + ...
            contract(Hpsi.identityChain, '12', contract(pairOpPsi, '23', pairOpPsi, '23*'), '13');
        % Add to opSumUp ID(1 ... i-1) x H(pair) and the same for
        % H.opSumDown
        Hpsi.opSumUp.next = contract(Hpsi.opSumUp.curr, '12', identityPair, '13');
        Hpsi.opSumDown.next = contract(Hpsi.opSumDown.curr, '12', identityPair, '13');
        Hpsi.opSumUp.next  = Hpsi.opSumUp.next + ...
            contract(Hpsi.identityChain, '12', pairOpPsiUp, '13');
        Hpsi.opSumDown.next  = Hpsi.opSumDown.next + ...
            contract(Hpsi.identityChain, '12', pairOpPsiDown, '13');
    end
    H.identityChain
    H.opSumUp.curr
    H.opSumDown.next
    H.opSumBoth.curr
    E2 = 0;