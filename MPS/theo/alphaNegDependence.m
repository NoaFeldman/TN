function res = alphaNegDependence(alphas, params, fixed)
    sigma = params(1);
    w0 = abs(params(2));
    w1 = abs(params(3));
    w2 = abs(params(4));
    res = w0 .* fluxNegsT0(alphas, fixed) .* a(alphas, 0, sigma).^4 .* a(alphas, 0, sigma/2) + ...
          w1 .* fluxNegsT1346(alphas, fixed) .* a(alphas, 0, sigma).^2 .* a(alphas, 0, sigma/2).^2 .* ...
                                          (a(alphas, 2*pi, sigma).^2 + a(alphas, -2*pi, sigma).^2) + ...
          w2 .* fluxNegsT25(alphas, fixed) .* a(alphas, 0, sigma).^4 .* ...
                                        (a(alphas, pi, sigma/2).^2 + a(alphas, -pi, sigma/2).^2);
end

function res = a(alphas, mu, sigma)
    res = exp(-((alphas - mu)/sigma).^2);
end