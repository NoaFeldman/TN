function A = getCkCqMatrix(psi)
    % Creates A_ij = <psi|Ck^dagger Cq|psi>
    M = getCiCjMatrix(psi, length(psi));
    S = realSpaceToDualSpace(length(psi));
    A = S' * M * S;
 end    
