function X = getOnePointFunc(t, epsilon)
    L = 600;
    vF = 2;
    X = abs(pi .* sinh(pi/L * 2 * epsilon) ./ ...
        (2 .* sinh(pi / L .* (epsilon + i .* vF .* t)) .* ...
         sinh(pi / L .* (epsilon - i.*vF.*t))));
%     X = (t.^2 + epsilon^2) / (2*epsilon);
end