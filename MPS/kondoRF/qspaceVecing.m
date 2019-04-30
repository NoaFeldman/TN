function [vecedLiou, vecedRho, id, idIn] = qspaceVecing(Liou, rho)
    % rho:  _                  _
    %   ->-|_|->-       >> ->-|_|-<-
    %
    % idDouble:
    %    _______|_______
    %   |  ___________  |
    %   | |     _     | |
    %   |_| ->-|_|-<- |_|
    %
    idIn = getIdentity(rho, 2, '-0');
    vecedRho = contract(rho, 1, idIn, '1*');
    id = getIdentity(vecedRho, 1, vecedRho, 2, 'rho');
    id = chargeConservingIdentity(id);
    
    % Assuming Liou, rho in full space (explicit zero blocks where needed)
    vecedLiou = QSpace();
    IDL = getIdentity(rho, 1);
    IDR = getIdentity(rho, 2);
    vecedLiou = vecedLiou + -1i * contract(contract(...
        qTranspose(IDR, idIn), 1, contract(Liou.H, 2, id, 1), 2, [2 1 3]), '12', id, '12*');
    vecedLiou = vecedLiou - -1i * contract(contract(...
        qTranspose(Liou.H, idIn), 1, contract(IDL, 2, id, 1), 2, [2 1 3]), '12', id, '12*');
    vecedLiou = vecedLiou + contract(contract(...
        qTranspose(Liou.L', idIn), 1, contract(Liou.L, 2, id, 1), 2, [2 1 3]), '12', id, '12*');
    vecedLiou = vecedLiou - 1/2 * contract(contract(...
        qTranspose(Liou.LdagL, idIn), 1, contract(IDL, 2, id, 1), 2, [2 1 3]), '12', id, '12*');
    vecedLiou = vecedLiou - 1/2 * contract(contract(...
        qTranspose(IDR, idIn), 1, contract(Liou.LdagL, 2, id, 1), 2, [2 1 3]), '12', id, '12*');
    
end

function id = chargeConservingIdentity(id)
    inds = [];
    for i = 1:length(id.Q{1})
        if isequal(id.Q{1}(i, :), -id.Q{2}(i, :))
            inds = [inds i];
        end
    end
    id = cutQSpaceRows(id, inds);
end

function res = qTranspose(O, idIn)
    res = contract(contract(idIn, '1*', O, 1), 2, idIn, 1);
end
