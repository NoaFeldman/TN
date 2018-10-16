function r = chargeRenyi3(cicj, Q, u1, v1, v2)
    % Based on B6 in https://arxiv.org/pdf/1804.00632.pdf
    % Assuming A1 and A2 are adjacent.
    u2 = v1+1;
    Gplus = getG(1, cicj, u1, v1, v2);
    Gminus = getG(-1, cicj, u1, v1, v2);
    I = eye (v2 - u1 + 1);
    Gp2 = Gplus^2;
    Gp3 = Gp2 * Gplus;
    GpGm = Gplus * Gminus;
    GmGp = Gminus * Gplus;
    Gp2Gm = Gp2 * Gminus;
    GpGmGp = GpGm * Gplus;
    L2 = v2 - u2 + 1;
    alphas = -3.14:1e-2:3.14;
    r = 0;
    for alpha = alphas
        r = r + exp(1i * alpha * Q) / (2 * pi) * ...
            fluxRenyi3(I, Gplus, Gminus, Gp2, Gp3, GpGm, GmGp, Gp2Gm, GpGmGp, alpha, L2);
    end
end

function r = fluxRenyi3(I, Gplus, Gminus, Gp2, Gp3, GpGm, GmGp, Gp2Gm, GpGmGp, alpha, L2)

    Tr = (-1/2 * det((I - 3.*Gplus + 3.*Gp2 - Gp3) / 2^3 + exp(1i .* alpha) .* (I + 3.*Gplus + 3.* Gp2 + Gp3) / 2^3) + ...
        1 * det((I - 2.*Gplus + Gp2 - Gminus + 2.*GpGm - Gp2Gm) / 2^3 + exp(1i .* alpha) .* (I + 2.*Gplus + Gp2 + Gminus + 2.*GpGm + Gp2Gm) / 2^3) +  ...
        1/2 * det((I - 2.* Gplus - Gminus + GmGp + Gp2 + GpGm - GpGmGp) / 2^3 +  exp(1i .* alpha) .* (I + 2.* Gplus + Gminus + GmGp + Gp2 + GpGm + GpGmGp) / 2^3));        
   
    r = Tr * exp(-1i * alpha * L2);
end