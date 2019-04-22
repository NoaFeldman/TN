function [L, basis, eVals, eStates] = arnoldi(rho0, Liou, K, applyLiou, calcOverlap, dir, newStateArr)
    % Implementation of the Arnoldi method as described in section 2.2 of
    % https://arxiv.org/pdf/1510.08634.pdf.
    % L here is the L matrix from the paper, Liou is the Liouvillian, basis is
    % the rho_i basis.
    L = zeros(K+1);
    basis{1} = rho0;
    for j = 1:K+1
        basis{j+1} = applyLiou(basis{j}, Liou);
        for i = 1:j
            L(i, j) = calcOverlap(basis{i}, basis{j+1});
            basis{j+1} = basis{j+1} - basis{i} * L(i, j);
        end
        L(j+1, j) = sqrt(calcOverlap(basis{j+1}, basis{j+1}));
        if L(j+1, j) < 1e-10
            K = j-1;
            break
        end
        basis{j+1} = basis{j+1}/L(j+1, j);
        basis{j+1} = reOrthonormalize(basis, j+1, calcOverlap);
    end
    L = L(1:K+1, 1:K+1);
    if strcmp(dir, 'right')
        [eVecs, eVals] = eig(L);
        eVals = diag(eVals);
        eStates = newStateArr(K+1);
        for i = 1:K+1
            for j = 1:K+1
                eStates(i) = eStates(i) + eVecs(j, i)*basis{j};
            end
        end
    elseif strcmp(dir, 'left')
        [~, eVals, eVecsL] = eig(L);
        eVecsL = eVecsL';
        eVals = diag(eVals);
        eStates = newStateArr(K);
        for i = 1:K
            for j = 1:K
                eStates(i) = eStates(i) + eVecsL(j, i)*basis{j};
            end
        end
    end
end

% Avoid accumulating error by making sure again that the new vector is
% orthogonalized to the former and that it is exactly normalized.
function vec = reOrthonormalize(basis, j, calcOverlap);
    vec = basis{j};
    for i = 1:j-1
        vec = vec - basis{i} * calcOverlap(basis{i}, vec);
%         vec.v = vec.v - calcOverlap(basis(i), vec) * basis(i).v;
    end
    vec = vec / sqrt(calcOverlap(vec, vec));
%     vec.v = vec.v / sqrt(calcOverlap(vec, vec));
end