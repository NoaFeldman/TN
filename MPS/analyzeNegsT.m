clear cs
clear bs
load('exactNegsQuench3.mat');
n = 3;
% data = load('exactNegsQuench12.mat');
% n = 1;

l = 1000;
rn = data.(strcat('r', int2str(n)));

zeroIndex = find(data.alphas == 0);
alphaRegion = 1:length(data.alphas) - 1 + 1;

for i = 5:length(data.t)    
    [a, ~, ~, ~, yfit] = fitnonlin(data.alphas(alphaRegion), data.alphas(alphaRegion), ...
        abs(rn(i, alphaRegion)/rn(i, zeroIndex)), 0.01.*data.alphas(alphaRegion), ...
        0.01.*abs(rn(i, alphaRegion)/rn(i, zeroIndex)), 'renyiNegT', [5 2 0.3], [l/2 data.t(i) n K]);
    plot(data.alphas, abs(rn(i, :)/rn(i, zeroIndex))); 
    hold on;
    plot(data.alphas(alphaRegion), yfit);
    xlabel('$\alpha$', 'Interpreter', 'latex');
    ylabel(strcat('$R_', int2str(n), '(\alpha)/R_', int2str(n), '$'), 'Interpreter', 'latex');
    hold off
    pause(0.1);
end