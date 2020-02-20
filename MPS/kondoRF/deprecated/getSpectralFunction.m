function As = getSpectralFunction(omegas, Sminus, rhoSS, rAlphas, lAlphas, lambdaAlphas, getOmegaAlphas)
    As = zeros(length(omegas), 1);
    omegaAlphas = getOmegaAlphas(Sminus, rhoSS, rAlphas, lAlphas);
    for o = 1:length(omegas)
        zAlphas = real(omegaAlphas) + ...
            imag(omegaAlphas).*(omegas(o) + imag(lambdaAlphas))./real(lambdaAlphas);
        As(o) = -1/pi * sum(zAlphas .* real(lambdaAlphas) ./ ...
            ((omegas(o) + imag(lambdaAlphas)).^2 + real(lambdaAlphas).^2));
    end
end