function res = rhoSinTheta(l, t, epsilon)
    % from (82) https://arxiv.org/pdf/1501.00568.pdf
        res = sign(l) .* 1i .* min(abs(l), t) .* (1 + epsilon^2 ./ (2 .* (min(abs(l), t).^2 - max(abs(l), t).^2)));
end