function r = chargeRenyi3(cicj, Q, u1, v1, v2)
    % Based on B6 in https://arxiv.org/pdf/1804.00632.pdf
    % Assuming A1 and A2 are adjacent.
    [I, Gplus, Gminus, Gp2, Gp3, GpGm, Gp2Gm, L2] = getRenyiNegNecessities(cicj, u1, v1, v2);
    stepsize = 1e-2;
    alphas = -3.14:stepsize:3.14;
    r = 0;
    for alpha = alphas
        r = r + exp(-1i * alpha * Q) / (2 * pi / stepsize) * ...
            fluxRenyi3(I, Gplus, Gminus, Gp2, Gp3, GpGm, Gp2Gm, alpha, L2);
    end
end