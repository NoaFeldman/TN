function SA = getEntanglementEntropy(t, a, L)
    % fitfun for fitgui/fitnonlin.
    % based on Eq. (89) in https://arxiv.org/pdf/0905.4013.pdf
    % and Eq. (36) in https://arxiv.org/pdf/1105.4846.pdf
    c = 1;
    SA = - c / 6 .* log(getOnePointFunc(t, a(1), L)) + a(2);
%     SA = - c / 12 .* log(getOnePointFunc(t, a(1), L)) + a(2);
end