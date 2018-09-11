function SA = getEntanglementEntropy(t, a, fixed)
    % fitfun for fitgui/fitnonlin.
    % based on Eq. (89) in https://arxiv.org/pdf/0905.4013.pdf
    % and Eq. (36) in https://arxiv.org/pdf/1105.4846.pdf
    L = fixed(1);
    model = fixed(2);
    c = 1;
    pointFunc = 1;
    SA = pointFunc * c / 6 .* log(getScaledVariable(t, a(1), L, model)) + a(2);
end 