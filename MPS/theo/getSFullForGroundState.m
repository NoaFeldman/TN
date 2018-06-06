function sFull = getSFullForGroundState(L, a, dummy)
    sFull =  1/6 .* log(onePointFunc(L, L/2)) + a(1);
end

function X = onePointFunc(L, l)
    X = 2 .* L ./ pi .* sin(pi .* l ./ L);
end