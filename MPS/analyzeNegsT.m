clear nucs
clear ws
clear epss
clear r0
clear r1
clear r2
clear rt0
clear rt1
clear rt2
load('exactNegQuenchFlux10000_1000-1000_3.mat');
n = 3;
% data = load('exactNegsQuench12.mat');
% n = 2;
K = 1;

l = 1000;
rn = data.(strcat('r', int2str(n)));

zeroIndex = find(data.alphas == 0);
alphaRegion = 1:length(data.alphas) - 1 + 1;

for i = 1:6
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
    rt0(i) = chargeNegT(data.t(i), [epss(i) ws(i) nucs(i)], [n K 0 l]);
    rt1(i) = chargeNegT(data.t(i), [epss(i) ws(i) nucs(i)], [n K 1 l]);
    rt2(i) = chargeNegT(data.t(i), [epss(i) ws(i) nucs(i)], [n K 2 l]);
end