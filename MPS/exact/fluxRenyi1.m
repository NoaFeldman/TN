function r = fluxRenyi1(I, Gplus, Gminus, alpha, L2)
% digits(128);
    Tr = 1\sqrt(2) * (exp(-1i * pi / 4) * det((I - Gplus)/2 + exp(1i * alpha) .* (I + Gplus)/2) + ...
                      exp(1i * pi / 4) * det((I - Gminus)/2 + exp(1i * alpha) .* (I + Gminus)/2));
    r = Tr * exp(-1i * alpha * L2);
end
