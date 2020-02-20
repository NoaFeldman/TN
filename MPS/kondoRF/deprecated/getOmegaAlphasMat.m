function omegaAlphas = getOmegaAlphasMat(Sminus, rhoSS, rAlphas, lAlphas)
    omegaAlphas = zeros(length(rAlphas(1, :)), 1);
    n = sqrt(length(rAlphas));
    if n == 0
        return
    end
    for alpha = 1:length(rAlphas(1, :))
        omegaAlphas(alpha) = trace(Sminus * reshape(rAlphas(:, alpha), [n n])) * ...
                             trace(reshape(lAlphas(:, alpha), [n n]) * ...
                             (Sminus'*reshape(rhoSS, [n n]) - reshape(rhoSS, [n n])*Sminus'));
    end
end
    