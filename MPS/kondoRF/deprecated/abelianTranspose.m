function res = abelianTranspose(qspace)
    % Transposition of a QSpace in abelian symmetries only!
    res = qspace;
    res.data = cellfun(@transpose, res.data, 'UniformOutput', 0);
end