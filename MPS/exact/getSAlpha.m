function s = getSAlpha(n, alpha, f, avg)
    if (nargin == 3)
        % Are we taking average value (for f acheived for \Delta N_A vs N_A) of alpha or sticking to GS eq.(9)
        avg = 0;
    end
    s = ones(1, length(alpha));
    for l = 1 : length(f)
        if (avg)
            s = s .* (exp(1i .* alpha/2) * f(l)^n + exp(-1i .* alpha/2)*(1 - f(l))^n);
        else
            s = s .* (exp(1i .* alpha) * f(l)^n + (1 - f(l))^n);
        end
    end
end