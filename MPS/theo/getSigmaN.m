function sigma = getSigmaN(t, a, fixed)
    % ON HOLD - PROBABLY DEPRACATED AFTER UPGRADE OF GS EQ(6)
    L = fixed(1);
    model = fixed(2);
    
    % fitfun for fitgui/fitnonlin.
    % based on file:///home/noa/Downloads/EEt_calc.pdf
    % and Eq. (36) in https://arxiv.org/pdf/1105.4846.pdf
    sigma = log(getScaledVariable(t, a(1), L, model)) ./ (2*pi^2) + a(2);
%     sigma = - 1/2.*log(getOnePointFunc(t, a(1), L)) ./ (2*pi^2) + a(2);
end