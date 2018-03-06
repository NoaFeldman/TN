function s = getSNA(n, f, x) 
    % S_1(N_A),  based on https://arxiv.org/pdf/1711.09418.pdf
    NSteps = 1e5;
    stepSize = 2 * pi / NSteps;
    s = zeros(1, length(x));
    for alpha = -pi : stepSize : pi
        sAlpha = getSAlpha(n, alpha, f);
        s = s + sAlpha .* exp(complex(0, -alpha .* x)) / (2 * pi / stepSize);
    end
end

function s = getSAlpha(n, alpha, f)
    s = 1;
    for l = 1 : length(f)
        s = s * (exp(1i * alpha) * f(l)^n + (1 - f(l))^n);
    end
end