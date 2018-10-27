function res = rhoCosTheta(l, t, epsilon)
    % from (82) https://arxiv.org/pdf/1501.00568.pdf
    res = max(abs(l), t) .* (1 + epsilon^2 ./ (2.*(max(abs(l), t).^2 - min(abs(l), t).^2)));
end