function [AV, A0, H0] = getHO(epsE, Ueh, U, epsH, omegaL, Omega)
    [F, Z, S, IS] = getLocalSpace('FermionS', 'Acharge,Aspin', 'NC', 1);
    AV = IS.E;
    AV.info.itags = {'sv', 'av*'};
    AC = getIdentity(AV, 2, IS.E, 1, 'ac');
    AC.info.itags{2} = 'sc';
    op = QSpace();
    % \sum_sigma epsE n_{e, \sigma} on condactance level
    nUp = contract(F(1)', '23', F(1), '13');
    nDown = contract(F(2)', '23', F(2), '13');
    nE = nUp + nDown;
    op = op + epsE * nE;
    % \sum_sigma Ueh n_h n_{e, \sigma} on conductance level
    nH = contract(F(2), '23', F(2)', '13');
    op = op + Ueh *  nH;
    % U n_up N_down
    intOp = contract(nUp, 2, nDown, 1);
    op = op + U * intOp;
    % (epsH - omegaL) nH
    op = op + (epsH - omegaL) * nH;
    H0 = contract(op, '12', contract(AC, '1*', AC, 1), '13');
    % Omega a_v,down a^dagger_c,dowm + HC
    vOp = contract(contract(AV, 1, F(2), 2), 2, AV, '1*', [3 1 2]);
    cOp = contract(contract(AC, 2, F(2)', 2), 3, AC, '2*', [4 1 3 5 2]);
    rabiOp = Omega * contract(vOp, '123', cOp, '123');
    H0 = H0 + rabiOp + rabiOp';
%     A0 = permute(getIdentity(H0, 1, IS.E, 1, 'a0'), [1 3 2]);
%     A0.info.itags{3} = 's0';
    A0 = permute(AC, [1 3 2]);
    AV.info.itags = {'sv', 'L00*'};
end