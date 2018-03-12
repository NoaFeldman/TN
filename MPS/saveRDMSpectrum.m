function saveRDMSpectrum(fileName, psi, truncErr)
    if (nargin == 2)
        truncErr = 0;
    end
    % Save all RDM eigenvalues of psi by spin (sub system is half the
    % lattice).
    % psi is expected to be left canonical.
    workPsi = psi;
    % Sweep to mid chain 
    k = length(workPsi);
    opts = {'Nkeep', 1024};
    while(k > length(workPsi)/2)
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

%     % TODO - just a test, delete and change back to row 38.
%     l = workPsi(1);
%     l = contract(getIdentity(real(l), 1, real(l), 2), '12*', l, '12');
%     for ind = 2:(length(workPsi)/2)
%         l = contract(l, 2, workPsi(ind), 1);
%         l = contract(getIdentity(real(l), 1, real(l), 2), '12*', l, '12');
%     end
%     r = workPsi(length(workPsi));
%     r = contract(r, '23', getIdentity(real(r), 2, real(r), 3), '12*');
%     for ind = (length(workPsi) - 1):-1:(length(workPsi)/2 + 1)
%         r = contract(workPsi(ind), 3, r, 1);
%         r = contract(r, '23', getIdentity(real(r), 2, real(r), 3), '12*');
%     end
%     M = contract(l, 2, r, 1);
%     
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