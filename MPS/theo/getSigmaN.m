function sigma = getSigmaN(t, a, L)
    % fitfun for fitgui/fitnonlin.
    % based on file:///home/noa/Downloads/EEt_calc.pdf
    % and Eq. (36) in https://arxiv.org/pdf/1105.4846.pdf
    sigma = - log(getOnePointFunc(t, a(1), L)) ./ (2*pi^2) + a(2);
end