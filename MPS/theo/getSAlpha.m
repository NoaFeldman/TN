function s = getSAlpha(n, alpha, f)
    s = ones(1, length(alpha));
    for l = 1 : length(f)
        s = s .* (exp(1i .* alpha) * f(l)^n + (1 - f(l))^n);
    end
end