data = load('exactNegsT0.mat');
n = 1;
rn = data.(strcat('r', int2str(n)));

zeroIndex = find(data.alphas == 0);
alphaRegion = 100:length(data.alphas) - 100;

% expected: c = K / (2 n pi^2)
f = fittype('exp(-x^2 * c)', 'independent', 'x', 'dependent', 'y');
for i = 1:length(data.lRights)
    [fg, gof] = fit(data.alphas(alphaRegion).', abs(rn(i, alphaRegion)/rn(i, zeroIndex)).', f, 'StartPoint', [0.1]);
%     plot(fg, data.alphas, abs(rn(i, :)/rn(i, zeroIndex))); 
%     xlabel('$\alpha$', 'Interpreter', 'latex');
%     ylabel(strcat('$R_', int2str(n), '(\alpha)/R_', int2str(n), '$'), 'Interpreter', 'latex');
%     pause(0.1);
    cs(i) = fg.c;
end

scatter(log(data.lRights.^2 ./ (data.lRights + data.lLeft)), cs);
hl = lsline;
B = [(hl.YData(1) - hl.YData(2))/(hl.XData(1) - hl.XData(2)), (hl.YData(2)*hl.XData(1) - hl.YData(1)*hl.XData(2))/(hl.XData(1)-hl.XData(2))];
disp(strcat('K = ', num2str(B(1) * 2 * n * pi^2)));
