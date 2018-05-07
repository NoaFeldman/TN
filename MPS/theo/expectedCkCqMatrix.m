function At = expectedCkCqMatrix(N, A0, t)
    At = A0;
    for k = 1 : length(At)
        for q = k : length(At)
            At(k, q) = At(k, q) * ...
                exp(1j * t * (getEk(N, k) - getEk(N, q)));
            At(q, k) = At(q, k) * ...
                exp(1j * t * (getEk(N, q) - getEk(N, k)));
        end
    end
end

function E = getEk(N, k)
    % Returns Ek for the k'th allowed value in the dual space (k is an
    % integer).
    % We here use open boundary conditions.
    E = 2 * cos(getMomentum(k, N));
end