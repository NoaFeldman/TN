function X = getOnePointFunc(t, epsilon, L)
    vF = 2*sin(pi * (L/2 + 1) / (L+1));
    X = abs(pi .* sinh(pi/L * 2 * epsilon) ./ ...
        (2 .* sinh(pi / L .* (epsilon + 1i .* vF .* t)) .* ...
         sinh(pi / L .* (epsilon - 1i.*vF.*t))));
%     X = abs((cos(pi.*(vF.*t +1i*epsilon)/L).^2 - cos(pi/16)^2) ... % l = L/16
%         .* (cos(2*pi/16) .* cosh(2*pi*epsilon / L) - cos(2*pi.*vF.*t/L) ...
%         +2*abs(cos(pi.*(vF.*t+1i*epsilon)/L).^2 - cos(pi/16)^2))).^-1;    
end