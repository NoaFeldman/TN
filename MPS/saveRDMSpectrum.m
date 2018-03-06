function saveRDMSpectrum(fileName, psi)
    % Save all RDM eigenvalues of psi by spin (sub system is half the
    % lattice).
    % psi is expected to be left canonical.
    workPsi = psi;
    % Sweep to mid chain 
    k = length(workPsi);
    opts = {'Nkeep', 1024};
    while(k > length(workPsi)/2)
        M = contract(workPsi(k-1), 3, workPsi(k), 1);
        [workPsi, truncErr] = decomposeAndTruncate(M, k-1, workPsi, '<<', opts);
        if (truncErr > 0)
            disp(truncErr);
        end
        k = k - 1;
    end
    % We take M to represent easily S^z_A:
    %      ______                ______
    %  ---|______|---   >>>  ===|______|---
    %       |  |                     |
    % Q on the double line is exactly S^z_A.
    M = contract(getIdentity(real(M), 1, real(M), 2), '12*', M, '12');
    % We now force the spin value of A by projection, and SVD
    % separately.
    [minSZ, maxSZ] = getMinMaxSZ(M);
    spectrum = containers.Map();
    for m = minSZ : maxSZ
        projected = project(M, m);
        if (~isempty(projected))
            [~, ~, I] = orthoQS(projected, 1, '<<', opts{:});
            spectrum(num2str(m)) = I.svd.^2;
        end
    end
    save(fileName, 'spectrum');
end
    
function projected = project(M, m)
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

function [minSZ, maxSZ] = getMinMaxSZ(M)
    % Returns minimum and maximum *absolute value* of S^z_A for M.
    minSZ = abs(M.Q{1}(1));
    maxSZ = abs(M.Q{1}(1));
    for i = 1:length(M.Q{1})
        if (abs(M.Q{1}(i)) > maxSZ)
            maxSZ = abs(M.Q{1}(i));
        elseif (abs(M.Q{1}(i)) < minSZ)
            minSZ = abs(M.Q{1}(i));
        end
    end
end