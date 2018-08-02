function S = getEE(f, x, N)
    % S(N_A),  based on https://arxiv.org/pdf/1711.09418.pdf (derivation by
    % n of eq. 2, 9)
    stepSize = 2 * pi / (N + 1);
    S = zeros(1, length(x));
    for alpha = stepSize : stepSize : stepSize  * N
        sAlpha = getSAlpha(1, alpha, f);
        dnSAlpha = 0;
        for l = 1 : length(f)
            dnSAlpha = dnSAlpha - sAlpha * getLDerivative(l, alpha, f) ...
                / getLDonation(l, 1, alpha, f);
        end
        S = S + dnSAlpha .* exp(complex(0, -alpha .* x)) / (2 * pi / stepSize);
    end
end


function s = getSAlpha(n, alpha, f)
    s = 1;
    for l = 1 : length(f)
        s = s * getLDonation(l, n, alpha, f);
    end
end

function res = getLDonation(l, n, alpha, f)
    res = exp(1i * alpha) * f(l)^n + (1 - f(l))^n;
end

function res = getLDerivative(l, alpha, f)
    if (f(l) ~= 0 && f(l) ~= 1)
        res = (exp(1i * alpha) * f(l) * log(f(l)) + (1 - f(l))*log((1 - f(l))));
    else
        res = 0;
    end
end