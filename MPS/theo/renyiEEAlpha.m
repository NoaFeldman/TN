function r = renyiEEAlpha(x, a, fixed)
    r = exp(-(x).^2./(2*a(2))) +  a(1) .* (exp(-(x - 2*pi).^2./(2*a(2))) + exp(-(x + 2*pi).^2/(2*a(2))));
end