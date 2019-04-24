function [L, basis, eVals, eStatesR, eStatesL] = arnoldi(rho0, Liou, K, applyLiou, calcOverlap, newStateArr, breakOpt)
    % Implementation of the Arnoldi method as described in section 2.2 of
    % https://arxiv.org/pdf/1510.08634.pdf.
    % L here is the L matrix from the paper, Liou is the Liouvillian, basis is
    % the rho_i basis.
    if ~exist('breakOpt','var')
        breakOpt = 'bySubSpace';
    end
    L = zeros(K+1);
    basis{1} = rho0;
    for j = 1:K+1
        basis{j+1} = applyLiou(basis{j}, Liou);
        for i = 1:j
            L(i, j) = calcOverlap(basis{i}, basis{j+1});
            basis{j+1} = basis{j+1} - basis{i} * L(i, j);
        end
        L(j+1, j) = sqrt(calcOverlap(basis{j+1}, basis{j+1}));
        if strcmp(breakOpt, 'bySubSpace') && (L(j+1, j) < 1e-10)
            K = j-1;
            break
        end
        basis{j+1} = basis{j+1}/L(j+1, j);
        basis{j+1} = reOrthonormalize(basis, j+1, calcOverlap);
    end
    L = L(1:K+1, 1:K+1);
    [eVecsR, eVals, eVecsL] = eig(L);
    eVals = diag(eVals);
    eStatesR = newStateArr(K+1);
    eStatesL = newStateArr(K+1);
    for i = 1:K+1
        for j = 1:K+1
            eStatesR(i) = eStatesR(i) + eVecsR(j, i)*basis{j};
            eStatesL(i) = eStatesL(i) + eVecsL(j, i)*basis{j};
        end
    end
%     for i = 1:length(eStatesL)
%         eStatesL = eStatesL / calcOverlap(
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

function vecs = orthoEVecs(vecs, calcOverlap)
    for j = 1:length(vecs)
        for i = 1:j-1
            vecs(j) = vecs(j) - calcOverlap(vecs(i), vecs(j)) .* vecs(i);
        end
        vecs(j) = vecs(j) / sqrt(calcOverlap(vecs(j), vecs(j)));
    end
end