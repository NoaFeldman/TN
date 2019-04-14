function [M, E0] = lanczos(HR, HL, HR2, HL2, H, k, psi)
    [T, base] = getTridiagonal(HR, HL, HR2, HL2, H, k, psi);
    % Numeric error accumulating here and take psi out of normalizaztion.
    if length(psi(1).info.itags) == 3
        bc = 'open';
    else
        bc = 'periodic';
    end
    [V, E] = eig(T);
    E = diag(E);
    [E0, V0] = min(E);
    M = QSpace;
    for i = 1 : length(V)
        M = M + V(i, V0) * base(i);
    end
    if strcmp(bc, 'open')
        norm = getscalar(contract(M, '1234', M, '1234*'));
    else
        norm = getscalar(contract(M, '123456', M, '123456*'));
    end
    for j = 1:length(M.data)
        M.data{j} = M.data{j}/sqrt(norm);
    end
end
        