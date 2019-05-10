function As = getSpectralFunctionQV(omegas, omegaAlphas, lambdaAlphas)
    for o = 1:length(omegas)
        zAlphas = real(omegaAlphas) + ...
            imag(omegaAlphas) .* (omegas(o) + imag(lambdaAlphas))./real(lambdaAlphas);
        As(o) = -1/pi .* sum(zAlphas .* real(lambdaAlphas) ./ ...
                             ((omegas(o) + imag(lambdaAlphas)).^2 + (real(lambdaAlphas)).^2));
    end
end