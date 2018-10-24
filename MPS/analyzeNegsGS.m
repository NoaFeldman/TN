clear cs
clear bs
data = load('exactNegsT0.mat');
n = 3;
rn = data.(strcat('r', int2str(n)));

zeroIndex = find(data.alphas == 0);
alphaRegion = 1:length(data.alphas) - 1 + 1;

% expected: c = K / (2 n pi^2)
% f = fittype('exp(-x^2 * c)', 'independent', 'x', 'dependent', 'y');
f = fittype(strcat('exp(-x^2 * (a + b + c)) + ', ...
    '0.5 * (exp(-x^2 * c - (x + 2*pi)^2 * a - (x - 2*pi)^2 * b) + ', ...
    'exp(-x^2 * c - (x + 2*pi)^2 * b - (x - 2*pi)^2 * a) + ', ...
    'exp(-x^2 * b - (x + 2*pi)^2 * c - (x - 2*pi)^2 * a) + ', ...
    'exp(-x^2 * b - (x + 2*pi)^2 * a - (x - 2*pi)^2 * c) + ', ...
    'exp(-x^2 * a - (x + 2*pi)^2 * b - (x - 2*pi)^2 * c) + ', ...
    'exp(-x^2 * a - (x + 2*pi)^2 * c - (x - 2*pi)^2 * b))'), ...
    'independent', 'x', 'dependent', 'y');

for i = 1:length(data.lRights)
%     [fg, gof] = fit(data.alphas(alphaRegion).', abs(rn(i, alphaRegion)/rn(i, zeroIndex)).', f, 'StartPoint', [0.1]);
    [fg, gof] = fit(data.alphas(alphaRegion).', abs(rn(i, alphaRegion)/rn(i, zeroIndex)).', f, 'StartPoint', [0.1 0.1 -0.1]);
    plot(fg, data.alphas, abs(rn(i, :)/rn(i, zeroIndex))); 
    xlabel('$\alpha$', 'Interpreter', 'latex');
    ylabel(strcat('$R_', int2str(n), '(\alpha)/R_', int2str(n), '$'), 'Interpreter', 'latex');
    pause(0.1);
    as(i) = fg.a;
    bs(i) = fg.b;
    cs(i) = fg.c;
  
end

% scatter(log(data.lRights.^2 ./ (data.lRights + data.lLeft)), cs); 
scatter(log(data.lRights.^2 ./ (data.lRights + data.lLeft)), as + bs + cs);
hl = lsline;
B = [(hl.YData(1) - hl.YData(2))/(hl.XData(1) - hl.XData(2)), (hl.YData(2)*hl.XData(1) - hl.YData(1)*hl.XData(2))/(hl.XData(1)-hl.XData(2))];
disp(strcat('K = ', num2str(B(1) * 2 * n * pi^2)));
