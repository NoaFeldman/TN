function res = abelianConj(qspace)
    res = qspace;
    res.data = cellfun(@conj, res.data, 'UniformOutput', 0);
end