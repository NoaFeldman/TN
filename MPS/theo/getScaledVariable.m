function X = getScaledVariable(var, epsilon, L, option)
    switch option
        case 1 % After quench, mid chain
            t = var;
            vF = 2*sin(pi * (L/2 + 1) / (L+1));
            % Half system (finite) after quantum quench in the middle:
            X = (abs(pi .* sinh(pi/L * 2 * epsilon) ./ ...
                (2 .* sinh(pi / L .* (epsilon + 1i .* vF .* t)) .* ...
                 sinh(pi / L .* (epsilon - 1i.*vF.*t))))).^(-1);
        case 2 % Ground state
            ratio = var;
            X = L.*sin(pi.*ratio);
        case 3 % Ground state, infinite system
            ratio = var;
            X = ratio * L;
        case 4 % After quench, infinite chain
            t = var;
            X = (epsilon^2 + t.^2) / (epsilon / 2);
        otherwise
            disp('Model type nonexistant.');
    end
%     lRatio = 0.25;
%     X = abs((cos(pi.*(vF.*t +1i*epsilon)/L).^2 - cos(pi*lRatio)^2) ... 
%         .* (cos(2*pi*lRatio) .* cosh(2*pi*epsilon / L) - cos(2*pi.*vF.*t/L) ...
%         -2*abs(cos(pi.*(vF.*t+1i*epsilon)/L).^2 - cos(pi*lRatio)^2))).^-1;    
end