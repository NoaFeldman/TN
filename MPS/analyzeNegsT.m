clear all
% load('exactNegsQuench3.mat');
% n = 3;
data = load('exactNegsQuench12.mat');
n = 1;

rn = data.(strcat('r', int2str(n)));

zeroIndex = find(data.alphas == 0);
% alphaRegion = 8:length(data.alphas) - 8;
alphaRegion = 100:length(data.alphas) - 100;

f = fittype('exp(-x^2 * c)', 'independent', 'x', 'dependent', 'y');
for i = 1:length(data.t)
    [fg, gof] = fit(data.alphas(alphaRegion).', abs(rn(i, alphaRegion)/rn(i, zeroIndex)).', f, 'StartPoint', [0.1]);
%     plot(fg, data.alphas, abs(rn(i, :)/rn(i, zeroIndex))); 
%     xlabel('$\alpha$', 'Interpreter', 'latex');
%     ylabel(strcat('$R_', int2str(n), '(\alpha)/R_', int2str(n), '$'), 'Interpreter', 'latex');
%     pause(0.1);
    cs(i) = fg.c;
end

%                               epsilon, const   l/v_F  K  
plot(renyiFluxSigma2(data.t, [0.1,    1],    [500    1  n]), cs);
pause(1);
[~, maxInd] = max(cs);
scatter(renyiFluxSigma2(data.t(1:maxInd), [0.1,    1],    [500    1  n]), cs(1:maxInd));

hl = lsline;
B = [(hl.YData(1) - hl.YData(2))/(hl.XData(1) - hl.XData(2)), (hl.YData(2)*hl.XData(1) - hl.YData(1)*hl.XData(2))/(hl.XData(1)-hl.XData(2))];
disp(strcat('1(?) = ', num2str(B(1) * (2 * pi)^2)));