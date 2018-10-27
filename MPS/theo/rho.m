function r = rho(l, t, epsilon)
    % from (82) https://arxiv.org/pdf/1501.00568.pdf
    r = ((epsilon^2 + l^2 + t.^2).^2 - 4 * l^2 .* t.^2).^(1/4);
end