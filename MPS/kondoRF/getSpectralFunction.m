function As = getSpectralFunction(omegas, omegaAlphas, lambdaAlphas)
    % Based on Eq. 18 from https://arxiv.org/pdf/1811.03518.pdf
    As = zeros(1, length(omegas));
    for a = 1:length(lambdaAlphas)
        if imag(lambdaAlphas(a)) == 0 || ... 
                abs(real(lambdaAlphas(a)) / imag(lambdaAlphas(a))) < 1e-16 % Equivalent for real(lambdaAlphas(a)) = 0
            % Use the Sokhotski–Plemelj theorem
            [~, peakInd] = min(abs(omegas + imag(lambdaAlphas(a))));
            As(peakInd) = As(peakInd) + omegaAlphas(a);
        else
            As = As - 1./pi .* imag(omegaAlphas(a)./(omegas + imag(lambdaAlphas(a)) - 1i .* real(lambdaAlphas(a))));
        end
    end
end