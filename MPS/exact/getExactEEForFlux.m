function S = getExactEEForFlux(f, alphas)
    % Some tests to check the alpha-dependency of the exact results.
    % S(alpha),  based on https://arxiv.org/pdf/1711.09418.pdf 
    
    % Some bad usage of Matlab's vector work in here, if I decide to use
    % this one more often, change.
    S = zeros(1, length(alphas));
    for a = 1:length(alphas)
        alpha = alphas(a);
        s1Alpha = getSAlpha(1, alpha, f);
        dnS1Alpha = 0;
        for l = 1 : length(f)
            dnS1Alpha = dnS1Alpha - s1Alpha * getLDerivative(l, alpha, f) ...
                / getLDonation(l, 1, alpha, f);
        end
        S(a) = -dnS1Alpha;
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