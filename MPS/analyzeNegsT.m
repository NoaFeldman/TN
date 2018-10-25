clear cs
clear bs
load('exactNegsQuench3.mat');
n = 3;
% data = load('exactNegsQuench12.mat');
% n = 1;

l = 1000;
rn = data.(strcat('r', int2str(n)));

zeroIndex = find(data.alphas == 0);
alphaRegion = 4:length(data.alphas) - 4 + 1;

f = fittype(strcat('exp(-x^2 * (b + b + c)) + ', ...
    '2 * 1.286 * exp(-x^2 * b - (x + 2*pi)^2 * b - (x - 2*pi)^2 * c) + ', ...
    '2 * 1.286 * exp(-x^2 * c - (x + 2*pi)^2 * b - (x - 2*pi)^2 * b) + ', ...
    '2 * 1.286 * exp(-x^2 * b - (x + 2*pi)^2 * c - (x - 2*pi)^2 * b)'), ...
    'independent', 'x', 'dependent', 'y');
for i = 1:length(data.t)
%     l12 =  analyticallyContinuedExpr(-l/2, data.t(i), 0.1)^(2/n);
%     l23 =  analyticallyContinuedExpr(l/2, data.t(i), 0.1)^(2/n);
%     l13 =  analyticallyContinuedExpr(0, data.t(i), 0.1)^(2/n);
%     etas = sqrt(eta13(l, data.t(i), 0.1).^(-8 / n) ./ eta12(l, data.t(i), 0.1).^(-8 /n));
%     f = fittype(strcat('(nuc^3 * ', num2str(etas * l12^2 * l23^2 / l13), ')^(-((x)/(2*pi))^2) + ', ...
%         'w * ((nuc^2 *', num2str(l12^2), ')^(-((x       )/(2*pi))^2) * (nuc^2 * ', ...
%                        num2str(l23^2), ')^(-((x + 2*pi)/(2*pi))^2) * (nuc^(-1) * ', ...
%                        num2str(l13^-1), ')^(-((x - 2*pi)/(2*pi))^2) + ', ...
%              '(nuc^2 *', num2str(l12^2), ')^(-((x       )/(2*pi))^2) * (nuc^2 * ', ...
%                        num2str(l23^2), ')^(-((x - 2*pi)/(2*pi))^2) * (nuc^(-1) * ', ...
%                        num2str(l13^-1), ')^(-((x + 2*pi)/(2*pi))^2) + ', ... 
%              '(nuc^2 *', num2str(l12^2), ')^(-((x + 2*pi)/(2*pi))^2) * (nuc^2 * ', ...
%                        num2str(l23^2), ')^(-((x       )/(2*pi))^2) * (nuc^(-1) * ', ...
%                        num2str(l13^-1), ')^(-((x - 2*pi)/(2*pi))^2) + ', ...
%              '(nuc^2 *', num2str(l12^2), ')^(-((x + 2*pi)/(2*pi))^2) * (nuc^2 * ', ...
%                        num2str(l23^2), ')^(-((x - 2*pi)/(2*pi))^2) * (nuc^(-1) * ', ...
%                        num2str(l13^-1), ')^(-((x       )/(2*pi))^2) + ', ...  
%              '(nuc^2 *', num2str(l12^2), ')^(-((x - 2*pi)/(2*pi))^2) * (nuc^2 * ', ...
%                        num2str(l23^2), ')^(-((x       )/(2*pi))^2) * (nuc^(-1) * ', ...
%                        num2str(l13^-1), ')^(-((x + 2*pi)/(2*pi))^2) + ', ...
%              '(nuc^2 *', num2str(l12^2), ')^(-((x - 2*pi)/(2*pi))^2) * (nuc^2 * ', ...
%                        num2str(l23^2), ')^(-((x + 2*pi)/(2*pi))^2) * (nuc^(-1) * ', ...
%                        num2str(l13^-1), ')^(-((x       )/(2*pi))^2))'), ...  
%                    'independent', 'x', 'dependent', 'y');

    
    [fg, gof] = fit(data.alphas(alphaRegion).', abs(rn(i, alphaRegion)/rn(i, zeroIndex)).', f, 'StartPoint', [1 1]);
    plot(fg, data.alphas, abs(rn(i, :)/rn(i, zeroIndex))); 
    xlabel('$\alpha$', 'Interpreter', 'latex');
    ylabel(strcat('$R_', int2str(n), '(\alpha)/R_', int2str(n), '$'), 'Interpreter', 'latex');
    pause(0.1);
%     nucs(i) = fg.nuc;
%     ws(i) = fg.w;
    bs(i) = fg.b;
    cs(i) = fg.c;
end

%                           epsilon, const   l/v_F  
plot(renyiFluxSigma2(data.t, [0.1,    1],    [l/2  n]), 2.*bs + cs);
cftRegion = 9:30; % Best for n = 3 
% cftRegion = 2:16;% Best for n = 1,2
scatter(renyiFluxSigma2(data.t(cftRegion), [0.1,    1],    [l/2 n]), 2.*bs(cftRegion) + cs(cftRegion));

hl = lsline;
B = [(hl.YData(1) - hl.YData(2))/(hl.XData(1) - hl.XData(2)), (hl.YData(2)*hl.XData(1) - hl.YData(1)*hl.XData(2))/(hl.XData(1)-hl.XData(2))];
disp(strcat('K = ', num2str(B(1) * (2 * pi)^2)));