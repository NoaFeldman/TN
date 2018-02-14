function At = expectedCkCqMatrix(psi0, A0, t)
    At = A0;
    for k = 1 : length(At)
        for q = k : length(At)
            At(k, q) = At(k, q) * ...
                exp(j * t * (getEk(length(psi0), k) - getEk(length(psi0), q)));
            At(q, k) = At(q, k) * ...
                exp(j * t * (getEk(length(psi0), q) - getEk(length(psi0), k)));
        end
    end
end

function E = getEk(N, k)
    % Returns Ek for the k'th allowed value in the dual space (k is an
    % integer).
    % We here use open boundary conditions.
    E = 2 * cos(getMomentum(k, N));
end