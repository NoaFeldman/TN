function [l, r, I] = myOrthoQS(M, idx, dir, opts)
    if (length(M.Q) == 3)
        [r, l, I] = orthoQS(M, idx, dir, opts{:});
    else
        [l, r, I] = orthoQS(M, idx, dir, opts{:});
    end
    if (I.Nkeep == I.Ntot)
        return
    end
    originalNkeep = I.Nkeep;
    Nkeep = originalNkeep;
    while (Nkeep < I.Ntot && abs(I.svd(originalNkeep) - I.svd(Nkeep + 1)) / I.svd(originalNkeep) < 1e-2)
        Nkeep = Nkeep + 1;
    end
    if (Nkeep ~= originalNkeep)
        [l, r, I] = orthoQS(M, idx, dir, 'Nkeep', Nkeep);
    end
end