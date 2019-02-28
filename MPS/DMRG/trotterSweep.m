function [psi, truncErr] = trotterSweep(trotterGates, psi, opts)
%     % Nearest neighbours Hamiltonian
%     truncErr = 0;
%     for k = length(psi) - 1 : -1 : 1
%         [err, psi] = applyHPair(trotterGates, k, psi, '<<', opts);
%         if (err > truncErr)
%             truncErr = err;
%         end
%     end
%     for k = 1 : length(psi) - 1
%         [err, psi] = applyHPair(trotterGates, k, psi, '>>', opts);
%         if (err > truncErr)
%             truncErr = err;
%         end
%     end
% end

    % Next nearest-neighbour Hamiltonian
    truncErr = 0;
    M = contract(psi(length(psi)), 1, psi(length(psi)-1), 3, [3 4 1 2]);
    for k = length(psi) : -1 : 2
        M = contract(M, '23', trotterGates.nearest(k-1), '12', [1 3 4 2]);
        if k > 2
            M = contract(M, 1, psi(k - 2), 3, [4 5 1 2 3]);
            M = contract(M, '24', trotterGates.next(k-2), '12', [1 4 2 5 3]);
            [l, r, I] = myOrthoQS(M, [1 2 3], '<<', opts);
            if (I.svd2tr > truncErr)
                truncErr = I.svd2tr;
            end
            psi(k) = QSpace(r);
            psi(k).info.itags{1} = strcat(int2str(k-1), 'a', psi(k).info.itags{1});
            M = QSpace(l);
            M.info.itags{4} = strcat(int2str(k-1), 'a', M.info.itags{4});
        else
            [l, r, I] = myOrthoQS(M, [1 2], '<<', opts);
            if (I.svd2tr > truncErr)
                truncErr = I.svd2tr;
            end
            psi(k) = QSpace(r);
            psi(k).info.itags{1} = strcat(int2str(k-1), 'a', psi(k).info.itags{1});
            psi(k-1) = QSpace(l);
            psi(k-1).info.itags{3} = strcat(int2str(k-1), 'a', psi(k-1).info.itags{3});
        end
    end
    
    M = contract(psi(1), 3, psi(2), 1);
    for k = 1:length(psi) - 1
        M = contract(M, '23', trotterGates.nearest(k), '12', [1 3 4 2]);
        if k < length(psi) - 1
            M = contract(M, 4, psi(k + 2), 1);
            M = contract(M, '24', trotterGates.next(k), '12', [1 4 2 5 3]);
            [l, r, I] = myOrthoQS(M, [1 2], '>>', opts);
            if (I.svd2tr > truncErr)
                truncErr = I.svd2tr;
            end
            psi(k) = QSpace(l);
            psi(k).info.itags{3} = strcat(int2str(k), 'a', psi(k).info.itags{3});
            M = QSpace(r);
            M.info.itags{1} = strcat(int2str(k), 'a', M.info.itags{1});
        else
            [l, r, I] = myOrthoQS(M, [1 2], '>>', opts);
            if (I.svd2tr > truncErr)
                truncErr = I.svd2tr;
            end
            psi(k) = QSpace(l);
            psi(k).info.itags{3} = strcat(int2str(k), 'a', psi(k).info.itags{3});
            psi(k+1) = QSpace(r);
            psi(k+1).info.itags{1} = strcat(int2str(k), 'a', psi(k+1).info.itags{1});
        end
    end
%     For imaginary time propagation, renormalize state.
%     psi = normalize(psi);
end       

function normalized = normalize(psi)
    c = getOverlap(psi, psi);
    for i = 1 : length(psi(1).data)
        psi(1).data{i} = psi(1).data{i} / sqrt(c);
    end
    normalized = psi;
end