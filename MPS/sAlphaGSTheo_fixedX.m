function s = sAlphaGSTheo_fixedX(alpha, params, fixed)
    % Testing, depricate soon
    X = fixed(3);
    w1 = params(1);
    w2 = params(2);
    n = fixed(1);
    k = fixed(2);
    s = X.^(-k/n.*(alpha/(2*pi)).^2) + ...
        w1.*X.^(-k/n.*((alpha + 2*pi)/(2*pi)).^2) + w1.*X.^(-k/n.*((alpha - 2*pi)/(2*pi)).^2) + ...
        w2.*X.^(-k/n.*((alpha + 4*pi)/(2*pi)).^2) + w2.*X.^(-k/n.*((alpha - 4*pi)/(2*pi)).^2);
end