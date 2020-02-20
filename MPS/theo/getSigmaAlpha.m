function sigma = getSigmaAlpha(t, a, fixed)
    a = real(a);
    L = fixed(1);
    model = fixed(2);
    K = fixed(3);
    pointFunc = fixed(4);
    epsilon = abs(a(1));
    u = fixed(5);
    sigma = 2 * pi^2 ./ (pointFunc * K .* (log(getScaledVariable(t, epsilon, L, model, u) * real(a(2)))));
end