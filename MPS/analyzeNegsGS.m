clear cs
clear bs
data = load('exactNegsT0.mat');
n = 3;
rn = data.(strcat('r', int2str(n)));

zeroIndex = find(data.alphas == 0);
alphaRegion = 1:length(data.alphas) - 1 + 1;

% expected: c = K / (2 n pi^2)
% f = fittype('exp(-x^2 * c)', 'independent', 'x', 'dependent', 'y');
f = fittype(strcat('exp(-x^2 * (b + b + c)) + ', ...
    '0.5 * (exp(-x^2 * c - (x + 2*pi)^2 * b - (x - 2*pi)^2 * b) + ', ...
           'exp(-x^2 * c - (x + 2*pi)^2 * b - (x - 2*pi)^2 * b) + ', ...
           'exp(-x^2 * b - (x + 2*pi)^2 * c - (x - 2*pi)^2 * b) + ', ...
           'exp(-x^2 * b - (x + 2*pi)^2 * b - (x - 2*pi)^2 * c) + ', ...
           'exp(-x^2 * b - (x + 2*pi)^2 * b - (x - 2*pi)^2 * c) + ', ...
           'exp(-x^2 * b - (x + 2*pi)^2 * c - (x - 2*pi)^2 * b))'), ...
           'independent', 'x', 'dependent', 'y');

for i = 1:length(data.lRights)
    l1 = data.lRights(i);
    l2 = data.lLeft;
    f = fittype(strcat('(nuc^3 * ', num2str(l1^2 * l2^2 / (l1 + l2)), ')^(-2 / 3 *(x/(2*pi))^2) + ', ... % n = 3.....
    '1.5e-4 * ((nuc^2 * ', num2str(l1^2), ')^(-2 / 3 *((x)/(2*pi))^2) * (nuc^2 * ', ...
            num2str(l2^2), ')^(-2 / 3 *((x + 2*pi)/(2*pi))^2) * (nuc^(-1) * ', num2str(1/(l1+l2)), ')^(-2 / 3 *((x - 2*pi)/(2*pi))^2) + ' , ...    
         '(nuc^2 * ', num2str(l1^2), ')^(-2 / 3 *((x)/(2*pi))^2) * (nuc^2 * ', ...
            num2str(l2^2), ')^(-2 / 3 *((x - 2*pi)/(2*pi))^2) * (nuc^(-1) * ', num2str(1/(l1+l2)), ')^(-2 / 3 *((x + 2*pi)/(2*pi))^2) + ' , ...
         '(nuc^2 * ', num2str(l1^2), ')^(-2 / 3 *((x + 2*pi)/(2*pi))^2) * (nuc^2 * ', ...
            num2str(l2^2), ')^(-2 / 3 *((x - 2*pi)/(2*pi))^2) * (nuc^(-1) * ', num2str(1/(l1+l2)), ')^(-2 / 3 *((x)/(2*pi))^2) + ' , ...    
         '(nuc^2 * ', num2str(l1^2), ')^(-2 / 3 *((x + 2*pi)/(2*pi))^2) * (nuc^2 * ', ...
            num2str(l2^2), ')^(-2 / 3 *((x)/(2*pi))^2) * (nuc^(-1) * ', num2str(1/(l1+l2)), ')^(-2 / 3 *((x - 2*pi)/(2*pi))^2) + ' , ...
         '(nuc^2 * ', num2str(l1^2), ')^(-2 / 3 *((x - 2*pi)/(2*pi))^2) * (nuc^2 * ', ...
            num2str(l2^2), ')^(-2 / 3 *((x + 2*pi)/(2*pi))^2) * (nuc^(-1) * ', num2str(1/(l1+l2)), ')^(-2 / 3 *((x)/(2*pi))^2) + ' , ...    
         '(nuc^2 * ', num2str(l1^2), ')^(-2 / 3 *((x - 2*pi)/(2*pi))^2) * (nuc^2 * ', ...
            num2str(l2^2), ')^(-2 / 3 *((x)/(2*pi))^2) * (nuc^(-1) * ', num2str(1/(l1+l2)), ')^(-2 / 3 *((x + 2*pi)/(2*pi))^2))'), ...       
        'independent', 'x');
    [fg, gof] = fit(data.alphas(alphaRegion).', abs(rn(i, alphaRegion)/rn(i, zeroIndex)).', f, 'StartPoint', [10]);
% fg.nuc = 120;
    plot(fg, data.alphas, abs(rn(i, :)/rn(i, zeroIndex))); 
    xlabel('$\alpha$', 'Interpreter', 'latex');
    ylabel(strcat('$R_', int2str(n), '(\alpha)/R_', int2str(n), '$'), 'Interpreter', 'latex');
    pause(0.1);
%     as(i) = fg.a;
    nucs(i) = fg.nuc;
%     ws(i) = fg.w; 
end

% scatter(log(data.lRights.^2 ./ (data.lRights + data.lLeft)), cs); 
scatter(log(data.lRights.^2 ./ (data.lRights + data.lLeft)), bs + cs);
hl = lsline;
B = [(hl.YData(1) - hl.YData(2))/(hl.XData(1) - hl.XData(2)), (hl.YData(2)*hl.XData(1) - hl.YData(1)*hl.XData(2))/(hl.XData(1)-hl.XData(2))];
disp(strcat('K = ', num2str(B(1) * 2 * n * pi^2)));
