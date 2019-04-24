function omegaAlphas = getOmegaAlphasMat(Sminus, rhoSS, rAlphas, lAlphas)
    n = sqrt(length(rhoSS));
    SplusCommutator = kron(eye(n), Sminus') - kron(transpose(Sminus'), eye(n));
    SminusVeced = kron(eye(n), Sminus);
    omegaAlphas = zeros(length(rAlphas(1, :)), 1);
    for alpha = 1:length(rAlphas(1, :))
        omegaAlphas(alpha) = trace(Sminus * reshape(rAlphas(:, alpha), [n n])) * ...
                             trace(reshape(lAlphas(:, alpha), [n n]) * ...
                             (Sminus'*reshape(rhoSS, [n n]) - reshape(rhoSS, [n n])*Sminus'));
    end
end
    