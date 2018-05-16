function saveRDMSpectrum(fileName, psi, subsystemASize, truncErr)
    if (nargin == 3)
        truncErr = 0;
    end
    % Save all RDM eigenvalues of psi by spin (sub system is half the
    % lattice).
    % psi is expected to be left canonical.
    workPsi = psi;
    % Sweep to mid chain 
    k = length(workPsi);
    opts = {'Nkeep', 1024};
    while(k > subsystemASize)
        M = contract(workPsi(k-1), 3, workPsi(k), 1);
        [workPsi, Err] = decomposeAndTruncate(M, k-1, workPsi, '<<', opts);
        if (Err > 1e-16)
            disp(Err);
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
        projected = projectSZ(M, m);
        if (~isempty(projected))
            [~, ~, I] = orthoQS(projected, 1, '<<', opts{:});
            spectrum(num2str(m)) = I.svd.^2;
        end
    end
    save(fileName, 'spectrum', 'truncErr');
end

function [minSZ, maxSZ] = getMinMaxSZ(M)
    % Returns minimum and maximum value of S^z_A for M.
    minSZ = M.Q{1}(1);
    maxSZ = M.Q{1}(1);
    for i = 1:length(M.Q{1})
        if (M.Q{1}(i) > maxSZ)
            maxSZ = M.Q{1}(i);
        elseif (M.Q{1}(i) < minSZ)
            minSZ = M.Q{1}(i);
        end
    end
end