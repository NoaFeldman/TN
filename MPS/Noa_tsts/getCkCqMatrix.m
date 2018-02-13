function A = getCkCqMatrix(psi)
    % Creates A_ij = <psi|Ck^dagger Cq|psi>
    M = getCiCjMatrix(psi);
    S = realSpaceToDualSpace(length(psi));
    A = S' * M * S;
 end    

function M = getCiCjMatrix(psi)
    % Creates M_ij = <psi|Ci^dagger Cj|psi>
    M = zeros(length(psi), length(psi));
    for i = 1 : length(M)
        for j = i : length(M)
            M(i, j) = getCiCj(i, j, psi);
            M(j, i) = conj(M(i, j));
        end
    end
end