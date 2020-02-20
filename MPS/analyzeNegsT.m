% K = getK(0.5); 
% l = 64;  
% L = 256; 
% dmrgNegs;
% rn = rAlpha;
 
K = getK(0); 
l = 64;
L = 256;
data = load('exactNegQuenchFlux256_64-64_largeTime.mat');
rn = data.(strcat('r', int2str(n)));
t = data.t;

alphas = -3.14:0.01:3.14;
zeroIndex = find(alphas == 0);
alphaRegion = 1:length(alphas);

dirName = strcat('negsT_n', int2str(n), 'K', num2str(K), 'L', int2str(L));
mkdir(dirName)

for i = 1:length(t)
    [a, ~, ~, chisq, yfit] = fitnonlin(alphas(alphaRegion), alphas(alphaRegion), ...
        abs(rn(i, alphaRegion)/rn(i, zeroIndex)), 0.01.*ones(length(alphaRegion), 1), ...
         0.01.*abs(rn(i, alphaRegion)/rn(i, zeroIndex)), 'fluxNeg1stOrder', [1e-1 10], [l n K L t(i)], []);
%     plot(alphas, real(rn(i, :)/rn(i, zeroIndex))); 
%     hold on;
%     [a, ~, ~, chisq, yfit] = fitnonlin(alphas(alphaRegion), alphas(alphaRegion), ...
%         real(rn(i, alphaRegion)/rn(i, zeroIndex)), 0.01.*ones(length(alphaRegion), 1), ...
%          0.01.*real(rn(i, alphaRegion)/rn(i, zeroIndex)), 'fluxNegT', [1 60 0.35], [l n K L t(i)], []);
%     plot(alphas(alphaRegion), yfit);
%     xlabel('$\alpha$', 'Interpreter', 'latex');
%     ylabel(strcat('$R_', int2str(n), '(\alpha)/R_', int2str(n), '$'), 'Interpreter', 'latex');
%     title(strcat('$t = ', num2str(t(i)), '$'), 'Interpreter', 'latex');
%     saveas(gcf, strcat(dirName, '/t', int2str(t(i)), '.png'));
%     hold off
%     pause(0.1);
%     chisqs(i) = chisq;
%     epss(i) = a(1);
%     nucs(i) = a(2);
%     wpm1s(i) = a(3);
    rEven(i) = sum((exp(-0i.*alphas) + 2.* (exp(-2i.*alphas) + exp(-4i.*alphas) + exp(-6i.*alphas) + exp(-8i.*alphas) + exp(-10i.*alphas))).*rn(i, :)/rn(i, zeroIndex)./(2* pi) * 1e-2);
    rOdd(i) = sum((2.* (exp(-1i.*alphas) + exp(-3i.*alphas) + exp(-5i.*alphas) + exp(-7i.*alphas) + exp(-9i.*alphas))).*rn(i, :)/rn(i, zeroIndex)./(2* pi) * 1e-2);
    rQ0(i) = sum(exp(-0i.*alphas).*rn(i, :)./(rn(i, zeroIndex)*2* pi) * 1e-2);
    rQ1(i) = sum(exp(-1i.*alphas).*rn(i, :)./(rn(i, zeroIndex)*2* pi) * 1e-2);
    rQ2(i) = sum(exp(-2i.*alphas).*rn(i, :)./(rn(i, zeroIndex)*2* pi) * 1e-2);
    rQ3(i) = sum(exp(-3i.*alphas).*rn(i, :)./(rn(i, zeroIndex)*2* pi) * 1e-2);
    rQ4(i) = sum(exp(-4i.*alphas).*rn(i, :)./(rn(i, zeroIndex)*2* pi) * 1e-2);
    rQ5(i) = sum(exp(-5i.*alphas).*rn(i, :)./(rn(i, zeroIndex)*2* pi) * 1e-2);
end

save(strcat(dirName, '/fitParams'), 'nucs', 'wpm1s', 'rQ0', 'rQ1', 'rQ2', 'rQ3', 'rQ4', 'rQ5', 'n', 'l', 'L', 'K', 't', 'rn', 'alphas');