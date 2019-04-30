function As = getSpectralFunctionQV(omegas, Liou, rhoSS, Sminus)
    [vecedLiou, vecedRho, id, idIn] = qspaceVecing(Liou, rhoSS);
    [rAlphas, lDagAlphas, lambdaAlphas] = diagVecedLiou(vecedLiou, id);
    omegaAlphas = getOmegaAlphasQV(rAlphas, lAlphas, Sminus, rhoSS, id, idIn);
    Zs = real(omegaAlphas)' + imag(omegaAlphas)' / real(lambdaAlphas) .* ...
                              (omegas' + imag(lambdaAlphas));
    % ...
end