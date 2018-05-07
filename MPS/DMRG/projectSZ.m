function projected = projectSZ(M, m)
    % M is a two-site matrix:
    %      ______
    % ====|______|--
    %          |
    % Projects S^z_A = m on M
    projector = getIdentity(M, 1);
    for i = 1:length(projector.Q{1})
        if (projector.Q{1}(i) == m)
            for j = 1:length(projector.Q)
                projector.Q{j} = projector.Q{j}(i);
            end
            projector.data = projector.data(i); 
            projected = contract(projector, 2, M, 1);
            return;
        end
    end
    projected = QSpace;
end
