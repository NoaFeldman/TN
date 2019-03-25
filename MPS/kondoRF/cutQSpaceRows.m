function qf = cutQSpaceRows(qi, keepRows)
    qf = qi;
    for i = 1:length(qi.Q)
        qf.Q{i} = qf.Q{i}(keepRows, :);
    end
    qf.data = qf.data(keepRows);
end