function omegaAlphas = getOmegaAlphasQSpace(Sminus, rhoSS, rAlphas, lAlphas)
    omegaAlphas = zeros(length(rAlphas), 1);
    leftOp = Sminus' * rhoSS - rhoSS * Sminus';
    for alpha = 1:length(rAlphas)
        omegaAlphas(alpha) = trace(Sminus * rAlphas(alpha)) * ...
                             trace(lAlphas(alpha) * leftOp);
    end
end