function s = getSNA(n, f, x, N) 
    % S_n(N_A),  based on https://arxiv.org/pdf/1711.09418.pdf
    stepSize = 2 * pi / (N + 1);
    s = zeros(1, length(x));
    for alpha = stepSize : stepSize : stepSize  * N
        sAlpha = getSAlpha(n, alpha, f);
        s = s + sAlpha .* exp(complex(0, -alpha .* x)) / (2 * pi / stepSize);
    end
end