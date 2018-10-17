function r = fluxRenyi3(I, Gplus, Gminus, Gp2, Gp3, GpGm, Gp2Gm, alpha, L2)
    Tr = -1/2 * det((I - 3.*Gplus + 3.*Gp2 - Gp3) / 2^3 + exp(1i .* alpha) .* (I + 3.*Gplus + 3.* Gp2 + Gp3) / 2^3) + ...
        3/2 * det((I - 2.*Gplus + Gp2 - Gminus + 2.*GpGm - Gp2Gm) / 2^3 + exp(1i .* alpha) .* (I + 2.*Gplus + Gp2 + Gminus + 2.*GpGm + Gp2Gm) / 2^3);
    r = Tr * exp(-1i * alpha * L2);
end