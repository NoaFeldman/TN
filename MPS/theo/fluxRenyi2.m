function r = fluxRenyi2(I, Gplus, Gminus, alpha, L2)
    Tr = det((I - Gminus)*(I - Gplus)/2^2 + exp(1i * alpha) .* (I + Gminus)*(I + Gplus)/2^2);
    r = Tr * exp(-1i * alpha * L2);
end