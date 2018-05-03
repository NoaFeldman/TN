function s = stheo(t, NA, epsS, c1, sigma2, L)
    % S_A(N_A) = (c_n * A * B * c)'
    % c_n' |_{n = 1} = c1 (from sFull fit)
    % A = X(t, epsS)^(c(n - 1/n)/6)
    % B = sqrt(pi * n / 2 * K * log(X(t, epsP)) = sqrt(n / (2 * pi *
    % sigma2))
    % C = e^(- n * N_A^2 / (2 * sigma2))
    global c;
    c = 1;
    cn = 1;
    xS = getOnePointFunc(t, epsS, L);
   
    s = - (cn .* dnA(xS) .* B(sigma2) .* C(NA, sigma2) + ...
        cn .* A(xS) .* dnB(sigma2) .* C(NA, sigma2) + ...
        cn .* A(xS) .* B(sigma2) .* dnC(NA, sigma2) - ...
        c1 .* A(xS) .* B(sigma2) .* C(NA, sigma2));
end
    
function a = A(x)
    a = 1;
end

function da = dnA(x)
    global c;
    da = c/6.* log(x);
end

function b = B(sigma2)
    b = sqrt(1 ./ (2 .* pi .* sigma2));
end

function db = dnB(sigma2)
    db = 1/2 .* sqrt(1 ./ (2 .* pi .* sigma2));
end

function c = C(NA, sigma2)
    c = exp(- NA^2 ./ (2 .* sigma2));
end

function dc = dnC(NA, sigma2)
    dc = - NA^2 ./ (2 .* sigma2) .* C(NA, sigma2);
end

