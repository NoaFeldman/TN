function der = wDerivative(l, t, epsilon, L)
%     if l == 0
%         der = 1 ./ sqrt(epsilon^2 + t.^2);
%     else
%         der = abs(1 ./ rho(l, t, epsilon) .* sqrt((l + rhoCosTheta(l, t, epsilon)).^2 + (1i .* t + rhoSinTheta(l, t, epsilon)).^2));
%     end
    vF = 2;
    theta = 1i .* vF .* t;
    if l == 0
        der = pi * sinh(2 * pi / L * epsilon)./(2 * L .* (sinh(pi / L .*( epsilon + theta)) .* sinh(pi / L .*( epsilon - theta))));
    else
        der = ((cos(pi/L .* (t + 1i * epsilon)).^2 - cos(pi/L .* l).^2) .* ...
            (cosh(2*pi/L .* l).*cosh(2*pi/L .* epsilon) - cos(2*pi/L .* t) + ...
                2.*abs(cos(pi/L .* (t - 1i* epsilon)).^2 - cos(pi/L .* l).^2)))^(-1/4);
    end
end