function [w11t, w22t, w33t, w12, w12t, w23, w23t, w13, w13t] = wDiffs(l, t, epsilon)
    % Analytically continued |omega_i - omega_j|, based on analytical calcs
    % and WCR.
    w11t = abs(2/epsilon .* (l + rhoCosTheta(l, t, epsilon)));
    w22t = abs(2/epsilon .* sqrt(epsilon^2 + t.^2));
    w33t = abs(2/epsilon .* abs(-l + rhoCosTheta(-l, t, epsilon)));
    w12  = abs(1/epsilon .* sqrt((l + rhoCosTheta(l, t, epsilon) - sqrt(epsilon^2 + t.^2)).^2 + rhoSinTheta(l, t, epsilon).^2));
    w12t = abs(1/epsilon .* sqrt((l + rhoCosTheta(l, t, epsilon) + sqrt(epsilon^2 + t.^2)).^2 + rhoSinTheta(l, t, epsilon).^2));
    w23  = abs(1/epsilon .* sqrt((-l + rhoCosTheta(-l, t, epsilon) - sqrt(epsilon^2 + t.^2)).^2 + rhoSinTheta(-l, t, epsilon).^2));
    w23t = abs(1/epsilon .* sqrt((-l + rhoCosTheta(-l, t, epsilon) + sqrt(epsilon^2 + t.^2)).^2 + rhoSinTheta(-l, t, epsilon).^2));
    w13  = abs(1/epsilon .* sqrt((2*l + rhoCosTheta(l, t, epsilon) - rhoCosTheta(-l, t, epsilon)).^2 + (rhoSinTheta(l, t, epsilon) - rhoSinTheta(-l, t, epsilon)).^2));
    w13t  = abs(1/epsilon .* sqrt((2*l + rhoCosTheta(l, t, epsilon) + rhoCosTheta(-l, t, epsilon)).^2 + (rhoSinTheta(l, t, epsilon) + rhoSinTheta(-l, t, epsilon)).^2));
end