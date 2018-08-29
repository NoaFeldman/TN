function s = getSNA(n, f, x, L) 
    % S_n(N_A),  based on https://arxiv.org/pdf/1711.09418.pdf
    stepSize = 2 * pi / (L + 1);
    s = zeros(1, length(x));
    for alpha = stepSize : stepSize : stepSize  * (L + 1)
        sAlpha = getSAlpha(n, alpha, f);
        s = s + sAlpha .* exp(complex(0, -alpha .* x)) / (2 * pi / stepSize);
    end
end