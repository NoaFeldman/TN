function overlap = getOverlap(psi, l, dir, overlap)
    if (strcmp(dir, '>>'))
        if (l == 1)
            overlap = contract(psi(1), '12', psi(1), '12*');
            return;
        end
        overlap = contract(contract(overlap, 1, psi(l), 1), '12', psi(l), '12*');
        return;
    end
