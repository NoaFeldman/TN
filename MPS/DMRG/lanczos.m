function [M, E0] = lanczos(HR, HL, HR2, HL2, H, k, psi)
    [T, base] = getTridiagonal(HR, HL, HR2, HL2, H, k, psi);
    % Numeric error accumulating here and take psi out of normalizaztion.
    if length(psi(1).info.itags) == 3
        bc = 'open';
    else
        bc = 'periodic';
    end
    for i = 1:length(base)
        for j = i+1:length(base)
            if strcmp(bc, 'open')
                base(i) = base(i) - getscalar(contract(base(i), '1234', base(j), '1234*'))*base(j);
            else
                base(i) = base(i) - getscalar(contract(base(i), '123456', base(j), '123456*'))*base(j);
            end
        end
    end
    for i = 1:length(base)
        if strcmp(bc, 'open')
            norm = getscalar(contract(base(i), '1234', base(i), '1234*'));
        else
            norm = getscalar(contract(base(i), '123456', base(i), '123456*'));
        end
        for j = 1:length(base(i).data)
            base(i).data{j} = base(i).data{j}/sqrt(norm);
        end
    end
    [V, E] = eig(T);
    E = diag(E);
    [E0, V0] = min(E);
    M = QSpace;
    for i = 1 : length(V)
        M = M + V(i, V0) * base(i);
    end
end
        