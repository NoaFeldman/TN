function der = wDerivative(l, t, epsilon)
    if l == 0
        der = 1 ./ sqrt(epsilon^2 + t.^2);
    else
        der = abs(1 ./ rho(l, t, epsilon) .* sqrt((l + rhoCosTheta(l, t, epsilon)).^2 + (1i .* t + rhoSinTheta(l, t, epsilon)).^2));
    end
end