function projectedEs = getProjectedEs(E)
    % project E for DN = 0, 1, 2
    projectedEs = QSpace(3);
    for DN = 0:2
        projectedEs(DN + 1) = E;
        Q = E.Q;
        indices = find(-Q{1} + Q{3} == 2*DN);
        for i = 1:length(Q)
            Q{i} = Q{i}(indices);
        end
        projectedEs(DN + 1).Q = Q;
        projectedEs(DN + 1).data = projectedEs(DN + 1).data(indices);
    end
end