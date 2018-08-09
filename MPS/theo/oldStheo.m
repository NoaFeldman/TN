function s = oldStheo(X_scaled, sigma2, NA, pointFunc)
    % ON HOLD - PROBABLY DEPRACATED AFTER UPGRADE OF GS EQ(6)

    % S_A(N_A) = (c_n * A * B * C)'
    % c_n' |_{n = 1} = c1 (from sFull fit)
    % A = X^(pointFunc * c(n - 1/n)/12)
    % B = sqrt( n / 2 * K * pi * sigma2) = sqrt(n / (2 * pi *
    % sigma2))
    % C = e^(- n * N_A^2 / (2 * sigma2))
    %
    % Note: I take c_n = 1, and add a multiplicative constant to X_scaled
    % detrmined by sFull.

    s = - (dnA(X_scaled, pointFunc) .* B(sigma2) .* C(NA, sigma2) + ...
        A(X_scaled) .* dnB(sigma2) .* C(NA, sigma2) + ...
        A(X_scaled) .* B(sigma2) .* dnC(NA, sigma2));    
end
    
function a = A(x)
    a = 1;
end

function da = dnA(x, pointFunc)
    c = 1;
    da = c/(6/pointFunc).* log(x);
end

function b = B(sigma2)
    b = sqrt(1 ./ (2 .* pi .* sigma2));
end

function db = dnB(sigma2)
    db = 1/2 .* sqrt(1 ./ (2 .* pi .* sigma2));
end

function c = C(NA, sigma2)
    c = exp(- NA.^2 ./ (2 .* sigma2));
end

function dc = dnC(NA, sigma2)
    dc = - NA.^2 ./ (2 .* sigma2) .* C(NA, sigma2);
end

