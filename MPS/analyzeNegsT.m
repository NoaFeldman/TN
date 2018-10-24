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

f = fittype(strcat('exp(-x^2 * (b + b + c)) + ', ...
    '2 * 10 * exp(-x^2 * b - (x + 2*pi)^2 * b - (x - 2*pi)^2 * c) + ', ...
    '2 * 10 * exp(-x^2 * c - (x + 2*pi)^2 * b - (x - 2*pi)^2 * b) + ', ...
    '2 * 10 * exp(-x^2 * b - (x + 2*pi)^2 * c - (x - 2*pi)^2 * b)'), ...
    'independent', 'x', 'dependent', 'y');
for i = 1:length(data.t)
    [fg, gof] = fit(data.alphas(alphaRegion).', abs(rn(i, alphaRegion)/rn(i, zeroIndex)).', f, 'StartPoint', [0.1 -0.1]);
    plot(fg, data.alphas, abs(rn(i, :)/rn(i, zeroIndex))); 
    xlabel('$\alpha$', 'Interpreter', 'latex');
    ylabel(strcat('$R_', int2str(n), '(\alpha)/R_', int2str(n), '$'), 'Interpreter', 'latex');
    pause(0.1);
    bs(i) = fg.b;
    cs(i) = fg.c;
%     ws(i) = fg.w;

%     hold on
%     w = 1;
%     nuConst = 100;
%     l1eff = scaledVar12(data.t(i), 0.1, [l/2 n]) * nuConst;
%     l2eff = scaledVar23(data.t(i), 0.1, [l/2 n]) * nuConst;
%     l12eff = scaledVar13(data.t(i), 0.1, [l/2 n]) * nuConst^(-2);
%     plot(data.alphas, (l1eff * l2eff * l12eff).^(-(data.alphas/(2*pi)).^2) + ...
%                    w.*(l1eff.^(-((data.alphas)/(2*pi)).^2) .* l2eff.^(-((data.alphas + 2*pi)/(2*pi)).^2) .* l12eff.^(-((data.alphas - 2*pi)/(2*pi)).^2) + ...
%                        l1eff.^(-((data.alphas + 2*pi)/(2*pi)).^2) .* l2eff.^(-((data.alphas - 2*pi)/(2*pi)).^2) .* l12eff.^(-((data.alphas)/(2*pi)).^2) + ...
%                        l1eff.^(-((data.alphas - 2*pi)/(2*pi)).^2) .* l2eff.^(-((data.alphas)/(2*pi)).^2) .* l12eff.^(-((data.alphas + 2*pi)/(2*pi)).^2))); 
%     hold off
end

%                           epsilon, const   l/v_F  
plot(renyiFluxSigma2(data.t, [0.1,    1],    [l/2  n]), 2.*bs + cs);
cftRegion = 9:30; % Best for n = 3 
% cftRegion = 2:16;% Best for n = 1,2
scatter(renyiFluxSigma2(data.t(cftRegion), [0.1,    1],    [l/2 n]), 2.*bs(cftRegion) + cs(cftRegion));

hl = lsline;
B = [(hl.YData(1) - hl.YData(2))/(hl.XData(1) - hl.XData(2)), (hl.YData(2)*hl.XData(1) - hl.YData(1)*hl.XData(2))/(hl.XData(1)-hl.XData(2))];
disp(strcat('K = ', num2str(B(1) * (2 * pi)^2)));