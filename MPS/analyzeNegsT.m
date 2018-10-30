clear nucs
clear ws
clear epss
clear r0
clear r1
clear r2
clear rt0
clear rt1
clear rt2
data = load('exactNegQuenchFlux10000_1000-1000_3.mat');
n = 3;
% data = load('exactNegsQuench12.mat');
% n = 2;
K = 1;

l = 1000;
rn = data.(strcat('r', int2str(n)));

zeroIndex = find(data.alphas == 0);
alphaRegion = 2:length(data.alphas) - 2 + 1;

for i = 1:length(data.t)
    [a, ~, ~, chisq, yfit] = fitnonlin(data.alphas(alphaRegion), data.alphas(alphaRegion), ...
        abs(rn(i, alphaRegion)/rn(i, zeroIndex)), 0.01.*ones(length(alphaRegion), 1), ...
         0.01.*ones(length(alphaRegion), 1), 'renyiNegT', [5 9 0.3], [l/2 data.t(i) n K]);
    plot(data.alphas, abs(rn(i, :)/rn(i, zeroIndex))); 
    hold on;
    plot(data.alphas(alphaRegion), yfit);
    xlabel('$\alpha$', 'Interpreter', 'latex');
    ylabel(strcat('$R_', int2str(n), '(\alpha)/R_', int2str(n), '$'), 'Interpreter', 'latex');
    title(strcat('$t = ', num2str(data.t(i)), '$'), 'Interpreter', 'latex');
    hold off
    pause(0.1);
    chisqs(i) = chisq;
    epss(i) = a(1);
    nucs(i) = a(2);
    ws(i) = a(3);
    r0(i) = sum(rn(i, :)/rn(i, zeroIndex)./(2* pi) * 1e-2);
    r1(i) = sum(exp(-1i.*data.alphas).*rn(i, :)/rn(i, zeroIndex)./(2* pi) * 1e-2);
    r2(i) = sum(exp(-2i.*data.alphas).*rn(i, :)/rn(i, zeroIndex)./(2* pi) * 1e-2);
    [w11t, w22t, w33t, w12, w12t, w23, w23t, w13, w13t] = wDiffs(l/2, data.t(i), epss(i));
    dw1 = wDerivative(l/2, data.t(i), epss(i));
    dw2 = wDerivative(0, data.t(i), epss(i));
    dw3 = wDerivative(-l/2, data.t(i), epss(i));
    % Find 10 factor!
    rt0(i) = chargeNegT(n, 0, K, dw1, dw2, dw3, w11t, w22t, w33t, w12, w12t, w23, w23t, w13, w13t, ws(i), abs(nucs(i)));
    rt1(i) = chargeNegT(n, 1, K, dw1, dw2, dw3, w11t, w22t, w33t, w12, w12t, w23, w23t, w13, w13t, ws(i), abs(nucs(i)));
    rt2(i) = chargeNegT(n, 2, K, dw1, dw2, dw3, w11t, w22t, w33t, w12, w12t, w23, w23t, w13, w13t, ws(i), abs(nucs(i)));
end