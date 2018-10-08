function r = chargeRenyi3(cicj, Q, u1, v1, v2)
    alphas = -pi:1e-3:pi;
    r = 0;
    for alpha = alphas
        r = r + exp(1i * alpha * Q) * fluxRenyi3(cicj, alpha, u1, v1, v2) / (2 * pi);
    end
end