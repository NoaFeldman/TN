function X = getOnePointFunc(t, epsilon, L)
    vF = 2*sin(pi * (L/2 + 1) / (L+1));
    w = 1i.*vF.*t;
    Y = abs(dZ(w, epsilon, L) ./ analyticContinuationForImZ(w, epsilon, L));
    %hopefully, substitue X with Y once its problems are solved.s
    X = abs(pi .* sinh(pi/L * 2 * epsilon) ./ ...
        (2 .* sinh(pi / L .* (epsilon + 1i .* vF .* t)) .* ...
         sinh(pi / L .* (epsilon - 1i.*vF.*t))));
end

% The following transformation is taken from Eq (34) in SD
% https://arxiv.org/pdf/1105.4846.pdf.
function zeta = zeta(w, epsilon, L)
    zeta = sqrt(sinh(pi/L .* (w + epsilon))./ sinh(pi/L .* (w - epsilon)));
end

function dZeta = dZeta(w, epsilon, L)
    numerator = sinh(pi/L .* (w + epsilon));
    dNumerator = pi/L .* cosh(pi/L ./ (w + epsilon));
    denumenator = sinh(pi/L .* (w - epsilon));
    dDenumenator = pi/L .* cosh(pi/L ./ (w - epsilon));
    dSqrtArg = (dNumerator.*denumenator - dDenumenator.*numerator)./ denumenator.^2;
    dZeta = 1/2 .* dSqrtArg.^(-1/2);
end

function z = z(w, epsilon, L)
    z = coth(pi*epsilon/(2*L)) .* (1 + zeta(w, epsilon, L)) ./ (1 - zeta(w, epsilon, L));
end

function dZ = dZ(w, epsilon, L)
    dZ = coth(pi*epsilon/(2*L)) .* 2 ./ (1 - zeta(w, epsilon, L)).^2 .* dZeta(w, epsilon, L);
end

% We start from real w and in the end take the analytical continuation to
% imaginary w. Thus Im(z) != imag(z), but the amount that would have been
% imaginary if w was real.
function im = analyticContinuationForImZ(w, epsilon, L)
    zeta1 = sqrt( - sinh(pi/L .* (w + epsilon))./ sinh(pi/L .* (w - epsilon)));
    im = coth(pi*epsilon/(2*L)) * 2.*zeta1./(1 + zeta1.^2);
end