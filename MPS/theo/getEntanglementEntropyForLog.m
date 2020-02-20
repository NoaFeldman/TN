function SA = getEntanglementEntropyForLog(logLeff, a, fixed)
    % fitfun for fitgui/fitnonlin.
    % based on Eq. (89) in https://arxiv.org/pdf/0905.4013.pdf
    % and Eq. (36) in https://arxiv.org/pdf/1105.4846.pdf
    L = fixed(1);
    model = fixed(2);
    c = 1;
    pointFunc = fixed(3);
    epsilon = 1e-4; % a(1);
    SA = pointFunc * c / 6 .* logLeff + a(1);
end 