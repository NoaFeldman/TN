function res = rhoCosTheta(l, t, epsilon)
    % from (82) https://arxiv.org/pdf/1501.00568.pdf
%     res = max(abs(l), t) .* (1 + epsilon^2 ./ (2.*(max(abs(l), t).^2 - min(abs(l), t).^2)));
    multSign = ones(1, length(t));
    if l > 0
        multSign(find(t.^2 > l^2 + epsilon^2)) = multSign(find(t.^2 > l^2 + epsilon^2)).*(-1);
    else
        multSign(find(t.^2 < l^2 + epsilon^2)) = multSign(find(t.^2 < l^2 + epsilon^2)).*(-1);
    end
    res = -rho(l, t, epsilon).*sqrt(1/2  + multSign.*(epsilon^2 + l.^2 + t.^2)./rho(l, t, epsilon).^2);
end