clear nucs
clear ws
data = load('exactNegsT0.mat');
n = 3;
K = 1;
rn = data.(strcat('r', int2str(n)));

zeroIndex = find(data.alphas == 0);
alphaRegion = 1:length(data.alphas) - 1 + 1;

% expected: c = K / (2 n pi^2)
% f = fittype('exp(-x^2 * c)', 'independent', 'x', 'dependent', 'y');

for i = 1:length(data.lRights)
    l1 = data.lRights(i);
    l2 = data.lLeft;
    [a, ~, ~, ~, yfit] = fitnonlin(data.alphas(alphaRegion), data.alphas(alphaRegion), ...
        abs(rn(i, alphaRegion)/rn(i, zeroIndex)), 0.01.*data.alphas(alphaRegion), ...
        0.01.*abs(rn(i, alphaRegion)/rn(i, zeroIndex)), 'renyiNegGS', [9 1], [l1 l2 n K]);
    plot(data.alphas, abs(rn(i, :)/rn(i, zeroIndex)));
    hold on
    plot(data.alphas(alphaRegion), renyiNegGS(data.alphas(alphaRegion), a, [l1 l2 n K]));    
    xlabel('$\alpha$', 'Interpreter', 'latex');
    ylabel(strcat('$R_', int2str(n), '(\alpha)/R_', int2str(n), '$'), 'Interpreter', 'latex');
    pause(0.1);
    hold off
    nucs(i) = a(1);
    ws(i) = a(2);
end
