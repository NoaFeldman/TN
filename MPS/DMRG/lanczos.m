function [M, E0] = lanczos(HR, HL, HR2, HL2, H, k, psi)
    %TODO verify normalized.
    [T, base] = getTridiagonal(HR, HL, HR2, HL2, H, k, psi);
    [V, E] = eig(T);
    E0 = E(1, 1);
    V0 = 1;
    for i = 2 : length(E)
        if (E(i, i) < E0)
            E0 = E(i, i);
            V0 = i;
        end
    end
    M = QSpace;
    for i = 1 : length(V)
        M = M + V(i, V0) * base(i);
    end
end
        