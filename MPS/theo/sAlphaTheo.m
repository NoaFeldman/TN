function s = sAlphaTheo(alpha, params, fixed)
    epsilon = params(1);
    multiplyingConst = params(2);
    w = fixed(4); % params(3);
    n = fixed(1);
    L = fixed(2);
    t = fixed(3);
    model = fixed(4);
    X = multiplyingConst .* getScaledVariable(t, epsilon, L, model);
    s = X.^(-1/n.*(alpha/(2*pi)).^2) + ...
        w.*X.^(-1/n.*((alpha + 2*pi)/(2*pi)).^2) + w.*X.^(-1/n.*((alpha - 2*pi)/(2*pi)).^2);
end