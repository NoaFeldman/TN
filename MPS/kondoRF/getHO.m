function [AV, A0, H0] = getHO(epsE, Ueh, U, epsH, omegaL, Omega)
    [F, Z, S, IS] = getLocalSpace('FermionS', 'Acharge,Aspin', 'NC', 1);
    
    % Get NRG site for the valence level
    AV = IS.E;
    AV.info.itags = {'av*', 'sv'};
    % Force up spin on valence level to be full (This is a redundent degree
    % of freedom for this problem)
    AV = cutQSpaceRows(AV, 3:4);
    
    % Get NRG site for conductance level
    AC = getIdentity(AV, 1, IS.E, 1, 'a0');
    AC.info.itags{2} = 'sc';
    AC = permute(AC, [1 3 2]);
    
    % Define QD Hamiltonian based on Eq. (1) at Sbierski et al
    HBase = contract(getIdentity(AV, 2, AC, 3), '3', ...
        getIdentity(AV, 2, AC, 3), '3*');

    H = QSpace();
    
    nH = contract(F(1), '13', F(1)', '23', [2 1]);
    nH.Q{1} = nH.Q{1}(1, :);
    nH.Q{2} = nH.Q{2}(1, :);
    nH.data = nH.data(1);
    nH.info.itags = {'sv', 'sv*'};
    H = H + (epsH - omegaL) * contract(nH, 2, HBase, 1);
       
    nUp = contract(F(1)', '23', F(1), '13');
    nDown = contract(F(2)', '23', F(2), '13');
    nE = nUp + nDown;
    nE.info.itags = {'sc', 'sc*'};
    H = H - Ueh * contract(contract(nH, 2, HBase, 1), 2, nE, 2, [1 4 2 3]);

    HC = epsE * nE + U * contract(nUp, 2, nDown, 1);
    H = H + contract(HBase, 4, HC, 1);
    
    eDown = F(2);
    eDown.info.itags = {'sc', 'sc*', 'm*'};
    hUp = F(2)';
    hUp.info.itags = {'sv', 'sv*', 'm'};
    hUpDagger = F(2);
    hUpDagger.info.itags = {'sv', 'sv*', 'm*'};
    rabiOp = contract(eDown, 3, hUp, 3, [3 1 4 2]);
    rabiOpDagger = contract(eDown', 3, hUpDagger, 3, [3 1 4 2]);
    H = H + Omega * (rabiOp + rabiOpDagger);
    
    % H0:
    %  ____    ____    ____
    % |_AV_|--|H_QD|--|_AV_|
    %  __|_   |    |   __|_
    % |_AC_|--|____|--|_AC_|
    %    |               |
    %
    % A0 = AC
    H0V = contract(AV, '2*', contract(AV, 2, H, 3), '2');
    H0 = contract(contract(AC, '13', H0V, '24'), '23', AC, '13*', [2 1]);
    A0 = AC;   
end