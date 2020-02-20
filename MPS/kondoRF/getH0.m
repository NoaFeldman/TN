function [AV, A0, H0, rabiOp] = getH0(epsE, Ueh, U, omegaDiff, Omega)
    [F, Z, S, IS] = getLocalSpace('FermionS', 'Acharge,Aspin', 'NC', 1);
    
    % Get NRG site for the valence level
    vac = getVac(IS.E);
    AV = getIdentity(vac, 1, IS.E, 1);
    AV.info.itags = {'ad', 'av*', 'sv'};
    % Force up spin on valence level to be full (This is a redundent degree
    % of freedom for this problem)
    AV = cutQSpaceRows(AV, 3:4);
    
    % Get NRG site for conductance level
    AC = getIdentity(AV, 2, IS.E, 1, 'a0');
    AC.info.itags{2} = 'sc';
    AC = permute(AC, [1 3 2]);
    
    % Define QD Hamiltonian based on Eq. (1) at Sbierski et al
    HBase = contract(getIdentity(AV, 3, AC, 3), '3', ...
        getIdentity(AV, 3, AC, 3), '3*');

    H = QSpace();
    
    nH = contract(F(1), '13', F(1)', '23', [2 1]);
    nH.Q{1} = nH.Q{1}(1, :);
    nH.Q{2} = nH.Q{2}(1, :);
    nH.data = nH.data(1);
    nH.info.itags = {'sv', 'sv*'};
    H = H + omegaDiff * contract(nH, 2, HBase, 1);
       
    nUp = contract(F(1)', '23', F(1), '13');
    nDown = contract(F(2)', '23', F(2), '13');
    nE = nUp + nDown;
    nE.info.itags = {'sc', 'sc*'};
    H = H - Ueh * contract(contract(nH, 2, HBase, 1), 2, nE, 2, [1 4 2 3]);

    HC = epsE * nE + U * contract(nUp, 2, nDown, 1);
    H = H + contract(HBase, 4, HC, 1);
    
    cDown = F(2);
    cDown.info.itags = {'sc', 'sc*', 'm*'};
    vDownDagger = F(2)';
    vDownDagger.info.itags = {'sv', 'sv*', 'm'};
    rabiOp = contract(vDownDagger, '23', ...
        contract(Z, 2, contract(cDown, 2, HBase, 2, [3 1 4 5 2]), 2, [2 1 3 4 5]), '15');
    H = H + Omega * (rabiOp + rabiOp');
    
    % H0:
    %  ____    ____    ____
    % |_AV_|--| H  |--|_AV_|
    %  __|_   |    |   __|_
    % |_AC_|--|____|--|_AC_|
    %    |               |
    %
    % A0 = AC
    H0V = contract(AV, '13*', contract(AV, 3, H, 3), '13');
    H0 = contract(contract(AC, '13', H0V, '24'), '23', AC, '13*', [2 1]);
    A0 = AC;
    H0 = H0 + 1e-99*getIdentity(A0, 2);
end

function vac = getVac(I)
    vac = I;
    vac = cutQSpaceRows(vac, 1);
    vac.Q{1} = [0 0];
    vac.Q{2} = [0 0];
end
