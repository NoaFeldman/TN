% clear nucs
% clear ws
% clear epss
clear r0
clear r1
clear r2
clear r3
% data = load('exactNegQuenchFlux10000_1000-1000_.mat');
n = 3;
K = 1;

l = 1000;
% rn = data.(strcat('r', int2str(n)));
load('rn3.mat');

zeroIndex = find(data.alphas == 0);
alphaRegion = 200:length(data.alphas) - 200 + 1;

for i = 1:length(t)
%     [a, ~, ~, chisq, yfit] = fitnonlin(data.alphas(alphaRegion), data.alphas(alphaRegion), ...
%         abs(rn(i, alphaRegion)/rn(i, zeroIndex)), 0.01.*ones(length(alphaRegion), 1), ...
%          0.01.*abs(rn(i, alphaRegion)/rn(i, zeroIndex)), 'fluxNegT', [17 1.1 0.1 0.55], [l n K t(i)], []);

    [a, ~, ~, chisq, yfit] = fitnonlin(data.alphas(alphaRegion), data.alphas(alphaRegion), ...
        abs(rn(i, alphaRegion)/rn(i, zeroIndex)), 0.01.*ones(length(alphaRegion), 1), ...
         0.01.*abs(rn(i, alphaRegion)/rn(i, zeroIndex)), 'fluxNeg1stOrder', [17 1.1 0.1], [l n K t(i)], []);
    plot(data.alphas, abs(rn(i, :)/rn(i, zeroIndex))); 
    hold on;
    plot(data.alphas(alphaRegion), yfit);
    [a, ~, ~, chisq, yfit] = fitnonlin(data.alphas(:), data.alphas(:), ...
        abs(rn(i, :)/rn(i, zeroIndex)), 0.01.*ones(length(data.alphas), 1), ...
         0.01.*abs(rn(i, :)/rn(i, zeroIndex)), 'fluxNegT', [17 a(2)*3.8 a(3)*3.8 4], [l n K t(i)], []);
    plot(data.alphas, yfit);
    xlabel('$\alpha$', 'Interpreter', 'latex');
    ylabel(strcat('$R_', int2str(n), '(\alpha)/R_', int2str(n), '$'), 'Interpreter', 'latex');
    title(strcat('$t = ', num2str(t(i)), '$'), 'Interpreter', 'latex');
    hold off
    pause(0.1);
    chisqs(i) = chisq;
    epss(i) = a(1);
    nucs(i) = a(2);
    nucs2(i) = a(3);
    ws(i) = a(4);
    r0(i) = sum(rn(i, :)/rn(i, zeroIndex)./(2* pi) * 1e-2);
    r1(i) = sum(exp(-1i.*data.alphas).*rn(i, :)/rn(i, zeroIndex)./(2* pi) * 1e-2);
    r2(i) = sum(exp(-2i.*data.alphas).*rn(i, :)/rn(i, zeroIndex)./(2* pi) * 1e-2);
    r3(i) = sum(exp(-3i.*data.alphas).*rn(i, :)/rn(i, zeroIndex)./(2* pi) * 1e-2);
    r4(i) = sum(exp(-4i.*data.alphas).*rn(i, :)/rn(i, zeroIndex)./(2* pi) * 1e-2);
    r5(i) = sum(exp(-5i.*data.alphas).*rn(i, :)/rn(i, zeroIndex)./(2* pi) * 1e-2);
end