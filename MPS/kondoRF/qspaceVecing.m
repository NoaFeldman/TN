function [vecedLiou, id, idIn] = qspaceVecing(Liou, rhoT)
    % rho:  _                  _
    %   ->-|_|->-    >>    ->-|_|-<-
    %
    % idDouble:
    %    _______|_______
    %   |  ___________  |
    %   | |     _     | |
    %   |_| ->-|_|-<- |_|
    %
    idIn = getIdentity(rhoT, 2, '-0');
    vecedRho = contract(rhoT, 1, idIn, '1*');
    id = getIdentity(vecedRho, 1, vecedRho, 2, 'rho');
    id = chargeConservingIdentity(id);
%     id = contract(id, 2, idIn, '1*');
    id = contract(id, 2, idIn', '1', [1 3 2]);
    
    % Assuming Liou, rho in full space (explicit zero blocks where needed)
    vecedLiou = QSpace();
    ID = getIdentity(rhoT, 1);
    vecedLiou = vecedLiou + -1i * contract(Liou.H, 1, contract(ID, 2, id, 1), 2, [2 1 3]);
    vecedLiou = vecedLiou - -1i * contract(ID, 1, contract(Liou.H, 2, id, 1), 2, [2 1 3]);
    vecedLiou = vecedLiou + contract(Liou.L', 1, contract(Liou.L, 2, id, 1), 2, [2 1 3]);
    vecedLiou = vecedLiou - 0.5 * contract(Liou.LdagL, 1, contract(ID, 2, id, 1), 2, [2 1 3]);
    vecedLiou = vecedLiou - 0.5 * contract(ID, 1, contract(Liou.LdagL, 2, id, 1), 2, [2 1 3]);
    
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

function res = qTranspose(O, idIn)
    res = contract(contract(idIn, '1*', O, 1), 2, idIn, 1);
end
