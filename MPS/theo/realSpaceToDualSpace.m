function S = realSpaceToDualSpace(N)
    % Creates the tranfer matrix from real space to dual space.
    % Cn = SnkCk
    % We assume open boundary conditions.
    S = zeros(N, N);
    for n = 1 : N
        for k = 1 : N
            S(n, k) = sin(n * getMomentum(k, N)) * sqrt(2 /(N+1));
        end
    end
end 