function res = analyticallyContinuedExpr(l, t, epsilon)
% abs(dw/dz / 2Re(w)), based on (24 - 25) in https://arxiv.org/pdf/1501.00568.pdf
    if l == 0
        res = epsilon ./ (2 .* (epsilon^2 + t.^2));
    else
        res = 1/2 .* 1./(rho(l, t, epsilon) .* (l + rhoCosTheta(l, t, epsilon))) .*...
            sqrt((l + rhoCosTheta(l, t, epsilon)).^2 + (1i .* t + rhoSinTheta(l, t, epsilon)).^2);        
    end
end

function r = rho(l, t, epsilon)
    r = ((epsilon^2 + l.^2 + t.^2).^2 - 4 * l.^2 .* t.^2).^(1/4);
end

function res = rhoCosTheta(l, t, epsilon)
    % from (82) https://arxiv.org/pdf/1501.00568.pdf
        res = max(abs(l), t) .* (1 + epsilon^2 ./ (2.*(max(abs(l), t).^2 - min(abs(l), t).^2)));
end

function res = rhoSinTheta(l, t, epsilon)
    % from (82) https://arxiv.org/pdf/1501.00568.pdf
        res = sign(l) .* 1i .* min(abs(l), t) .* (1 + epsilon^2 ./ (2 .* (min(abs(l), t).^2 - max(abs(l), t).^2)));
end