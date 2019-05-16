function [vecedLiou, id, idIn] = qspaceVecing(Liou, rhoT)
    %  rho  _         >>    vecedRho  _
    %   ->-|_|->-     >>          ->-|_|-<-
    %
    idIn = getIdentity(rhoT, 1, '-0');
    vecedRho = contract(idIn, '1*', rhoT, 1);
    %  id     |
    %      ___|____
    %     |        |
    %     |        |
    %     v        ^
    % 
    % ChargeConservingIdentity throws away blocks that don't conserve q1 +
    % q2 .
    id = getIdentity(vecedRho, 1, vecedRho, 2, 'rho');
    id = chargeConservingIdentity(id);
    id = contract(idIn', '1', id, 1);
    
    % Assuming Liou, rho in full space (explicit zero blocks where added)
    vecedLiou = QSpace();
    % [H, rho]:
    %  _             _
    % |H|           |H|
    %  |    | - |    |
    %  |_id_|   |_id_|
    %     |        |
    vecedLiou = vecedLiou + -1i * ...
        (contract(Liou.H, 1, id, 1) - contract(id, 2, Liou.H, 2, [1 3 2]));
    % L rho L':
    %  _    __
    % |L|  |L'|
    %  |    | 
    %  |_id_| 
    %     |   
    vecedLiou = vecedLiou + ...
        contract(contract(Liou.L, 1, id, 1), 2, Liou.L', 2, [1 3 2]);
    % {L'L, rho}:
    vecedLiou = vecedLiou - 0.5 * ...
        (contract(Liou.LdagL, 1, id, 1) + contract(id, 2, Liou.LdagL, 2, [1 3 2]));   
end

function id = chargeConservingIdentity(id)
    inds = [];
    for i = 1:size(id.Q{3}, 1)
        if isequal(id.Q{3}(i, :), [0 0])
            inds = [inds i];
        end
    end
    id = cutQSpaceRows(id, inds);
end