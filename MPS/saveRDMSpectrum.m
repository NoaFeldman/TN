function saveRDMSpectrum(N, Delta)
    % Save all RDM eigenvalues, separated by XXZ model Delta and by 
    % the half-system's spin.
    tic;

    opts = {'Nkeep', 1024, 'stol', 1e-5};
    [psi, H, HR, HL] = getGroundState(N, 0, 1, Delta, 0, opts);
    
    % Sweep to mid chain 
    i = length(psi);
    while(i > length(psi)/2)
        [HL, HR, psi, ~, i, ~] = dmrgStep(HL, HR, H, psi, i, '<<', opts);
    end
    M = contract(psi(length(psi)/2), 3, psi(length(psi)/2 + 1), 1);
    % We take M to represent easily S^z_A:
    %      ______                ______
    %  ---|______|---   >>>  ===|______|---
    %       |  |                     |
    % Q on the double line is exactly S^z_A.
    M = contract(getIdentity(M, 1, M, 2), '12*', M, '12');
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
    disp(strcat('Finished calculating GS and spectrum for lambda = ', num2str(Delta)));
    toc;
    save(strcat('spectrumN', int2str(N), 'D', abs(num2str(Delta)), '.mat'), 'spectrum');
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