function sigma = getSigmaAlpha(t, a, fixed)
    % ON HOLD - PROBABLY DEPRACATED AFTER UPGRADE OF GS EQ(6)
    L = fixed(1);
    model = fixed(2);
    K = fixed(3);
    pointFunc = fixed(4);
    epsilon = abs(a(1));
    sigma = 2 * pi^2 ./ (pointFunc * K .* (log(getScaledVariable(t, epsilon, L, model) * real(a(2)))));
end