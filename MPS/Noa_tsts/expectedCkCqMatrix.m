function At = expectedCkCqMatrix(psi, A, t)
    At = A;
    for k = 1 : length(At)
        for q = k : length(At)
            At(k, q) = At(k, q) * ...
                exp(j * t * (getEk(length(psi), q) - getEk(length(psi), k)));
            At(q, k) = At(q, k) * ...
                exp(j * t * (getEk(length(psi), q) - getEk(length(psi), k)));
        end
    end
end

function E = getEk(N, k)
    % Returns Ek for the k'th allowed value in the dual space (k is an
    % integer).
    % We here use open boundary conditions.
    E = 2 * cos(getMomentum(k, N));
end