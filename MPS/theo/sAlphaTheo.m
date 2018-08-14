function s = sAlphaTheo(alpha, params, fixed)
    epsilon = params(1);
    multiplyingConst = params(2);
    w1 = params(3);
    w2 = params(4);
    n = fixed(1);
    L = fixed(2);
    t = fixed(3);
    k = fixed(4);
    model = fixed(4);
    X = multiplyingConst .* getScaledVariable(t, epsilon, L, model);
    s = X.^(-k/n.*(alpha/(2*pi)).^2) + ...
        w1.*X.^(-k/n.*((alpha + 2*pi)/(2*pi)).^2) + w1.*X.^(-k/n.*((alpha - 2*pi)/(2*pi)).^2) + ...
        w2.*X.^(-k/n.*((alpha + 4*pi)/(2*pi)).^2) + w2.*X.^(-k/n.*((alpha - 4*pi)/(2*pi)).^2);
end