% clear nucs
% clear ws
% clear epss
data = load('exactNegQuench10000_1000_1000_full.mat');
n = 3;
K = getK(0); 

vF = 2;
l = 1000 / vF; % 32; %
L = 10000; % 256; %
% rn = rAlpha; 
rn = data.(strcat('r', int2str(n)));
t = data.t;

alphas = -3.14:0.01:3.14;
zeroIndex = find(alphas == 0);
alphaRegion = 30:length(alphas) - 30 + 1;

for i = 1:length(t)
%     [a, ~, ~, chisq, yfit] = fitnonlin(data.alphas(alphaRegion), data.alphas(alphaRegion), ...
%         abs(rn(i, alphaRegion)/rn(i, zeroIndex)), 0.01.*ones(length(alphaRegion), 1), ...
%          0.01.*abs(rn(i, alphaRegion)/rn(i, zeroIndex)), 'fluxNeg1stOrder', [17 1.1], [l n K t(i)], []);
    plot(alphas, abs(rn(i, :)/rn(i, zeroIndex))); 
    hold on;
%     plot(data.alphas(alphaRegion), yfit);
    [a, ~, ~, chisq, yfit] = fitnonlin(alphas(:), alphas(:), ...
        abs(rn(i, :)/rn(i, zeroIndex)), 0.01.*ones(length(alphas), 1), ...
         0.01.*abs(rn(i, :)/rn(i, zeroIndex)), 'fluxNegT', [50 20 0.9], [l n K L t(i)], []);
    plot(alphas, yfit);
    xlabel('$\alpha$', 'Interpreter', 'latex');
    ylabel(strcat('$R_', int2str(n), '(\alpha)/R_', int2str(n), '$'), 'Interpreter', 'latex');
    title(strcat('$t = ', num2str(t(i)), '$'), 'Interpreter', 'latex');
    hold off
    pause(0.1);
    chisqs(i) = chisq;
    epss(i) = a(1);
    nucs(i) = a(2);
    ws(i) = a(3);
    rEven(i) = sum((exp(-0i.*data.alphas) + 2.* (exp(-2i.*data.alphas) + exp(-4i.*data.alphas) + exp(-6i.*data.alphas) + exp(-8i.*data.alphas) + exp(-10i.*data.alphas))).*rn(i, :)/rn(i, zeroIndex)./(2* pi) * 1e-2);
    rOdd(i) = sum((2.* (exp(-1i.*data.alphas) + exp(-3i.*data.alphas) + exp(-5i.*data.alphas) + exp(-7i.*data.alphas) + exp(-9i.*data.alphas))).*rn(i, :)/rn(i, zeroIndex)./(2* pi) * 1e-2);
    rQ0(i) = sum(exp(-0i.*data.alphas).*rn(i, :)/rn(i, zeroIndex)./(2* pi) * 1e-2);
    rQ1(i) = sum(exp(-1i.*data.alphas).*rn(i, :)/rn(i, zeroIndex)./(2* pi) * 1e-2);
    rQ2(i) = sum(exp(-2i.*data.alphas).*rn(i, :)/rn(i, zeroIndex)./(2* pi) * 1e-2);
    rQ3(i) = sum(exp(-3i.*data.alphas).*rn(i, :)/rn(i, zeroIndex)./(2* pi) * 1e-2);
    rQ4(i) = sum(exp(-4i.*data.alphas).*rn(i, :)/rn(i, zeroIndex)./(2* pi) * 1e-2);
    rQ5(i) = sum(exp(-5i.*data.alphas).*rn(i, :)/rn(i, zeroIndex)./(2* pi) * 1e-2);
end