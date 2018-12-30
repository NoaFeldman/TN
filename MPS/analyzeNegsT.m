% clear nucs
% clear ws
% clear epss
% data = load('exactNegQuench10000_1000_1000_full.mat');
n = 1;
K = getK(-0.5); 

vF = 2;
l = 64; % 1000 / vF; % 
L = 256; % 10000; % 
rn = rAlpha; 
ts = [10 18 28 38 48 58];
% rn = data.(strcat('r', int2str(n)));
% t = data.t;

alphas = -3.14:0.01:3.14;
zeroIndex = find(alphas == 0);
alphaRegion = 200:length(alphas) - 200 + 1;

for i = 1:length(ts)
    [a, ~, ~, chisq, yfit] = fitnonlin(data.alphas(alphaRegion), data.alphas(alphaRegion), ...
        abs(rn(i, alphaRegion)/rn(i, zeroIndex)), 0.01.*ones(length(alphaRegion), 1), ...
         0.01.*abs(rn(i, alphaRegion)/rn(i, zeroIndex)), 'fluxNeg1stOrder', [1e-1 10], [l n K L t(i)], []);
    plot(alphas, real(rn(i, :)/rn(i, zeroIndex))); 
    hold on;
    plot(data.alphas(alphaRegion), yfit);
%     [a, ~, ~, chisq, yfit] = fitnonlin(alphas(:), alphas(:), ...
%         real(rn(i, :)/rn(i, zeroIndex)), 0.01.*ones(length(alphas), 1), ...
%          0.01.*real(rn(i, :)/rn(i, zeroIndex)), 'fluxNegT', [1 1 0.1], [l n K L ts(i)], []);
%     plot(alphas, yfit);
    xlabel('$\alpha$', 'Interpreter', 'latex');
    ylabel(strcat('$R_', int2str(n), '(\alpha)/R_', int2str(n), '$'), 'Interpreter', 'latex');
    title(strcat('$t = ', num2str(ts(i)), '$'), 'Interpreter', 'latex');
    hold off
    pause(0.1);
    chisqs(i) = chisq;
    epss(i) = a(1);
    nucs(i) = a(2);
%     ws(i) = a(3);
%     rEven(i) = sum((exp(-0i.*data.alphas) + 2.* (exp(-2i.*data.alphas) + exp(-4i.*data.alphas) + exp(-6i.*data.alphas) + exp(-8i.*data.alphas) + exp(-10i.*data.alphas))).*rn(i, :)/rn(i, zeroIndex)./(2* pi) * 1e-2);
%     rOdd(i) = sum((2.* (exp(-1i.*data.alphas) + exp(-3i.*data.alphas) + exp(-5i.*data.alphas) + exp(-7i.*data.alphas) + exp(-9i.*data.alphas))).*rn(i, :)/rn(i, zeroIndex)./(2* pi) * 1e-2);
%     rQ0(i) = sum(exp(-0i.*data.alphas).*rn(i, :)/rn(i, zeroIndex)./(2* pi) * 1e-2);
%     rQ1(i) = sum(exp(-1i.*data.alphas).*rn(i, :)/rn(i, zeroIndex)./(2* pi) * 1e-2);
%     rQ2(i) = sum(exp(-2i.*data.alphas).*rn(i, :)/rn(i, zeroIndex)./(2* pi) * 1e-2);
%     rQ3(i) = sum(exp(-3i.*data.alphas).*rn(i, :)/rn(i, zeroIndex)./(2* pi) * 1e-2);
%     rQ4(i) = sum(exp(-4i.*data.alphas).*rn(i, :)/rn(i, zeroIndex)./(2* pi) * 1e-2);
%     rQ5(i) = sum(exp(-5i.*data.alphas).*rn(i, :)/rn(i, zeroIndex)./(2* pi) * 1e-2);
end