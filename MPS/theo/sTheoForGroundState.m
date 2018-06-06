function s = sTheoForGroundState(L, l, NA, onePointFunc)
    % S_A(N_A) = (A * B * c)'
    % A = X(t, epsS)^(c(n - 1/n)/6)
    % B = sqrt(pi * n / 2 * K * log(X(t, epsP)) = sqrt(n / (2 * pi *
    % sigma2))
    % C = e^(- n * N_A^2 / (2 * sigma2))
    constP = 0.577 + 1 + log(2);
    sigma2 = (log(onePointFunc(L, L/2)) + constP)/(2*pi^2);
%     constS = 0.726/3;
    constS = 0.726*6;
    xS = exp(-constS) .* onePointFunc(L, L/2);
    s = (dnA(xS) .* B(sigma2) .* C(NA, sigma2) + ...
        A(xS) .* dnB(sigma2) .* C(NA, sigma2) + ...
        A(xS) .* B(sigma2) .* dnC(NA, sigma2));    
end
    
function a = A(x)
    a = 1;
end

function da = dnA(x)
    da = 1/6.* log(x);
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

    