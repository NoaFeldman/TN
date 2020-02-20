data = load('exactQuench10000_0-5000_.mat')
L = 10000;
l = 1000;
K = 1;
n = 3;
zeroIndex = find(data.alphas == 0);
f = fittype('w1 * exp(-(x)^2 * s) +  w2 * (exp(-(x - 2*pi)^2 * s) + exp(-(x + 2*pi)^2 * s))', 'independent', 'x', 'dependent', 'y');
% for i = 1:length(data.t)
%     [fg, gof] = fit(data.alphas.', real(data.s3Alpha(:, i)./data.s3Alpha(zeroIndex, i)), ...
%         f, 'StartPoint', [0.5 1 0.8]);
% %     if (mod(i, 10) == 0)
% %         plot(fg, data.alphas, real(data.s3Alpha(:, i)./data.s3Alpha(zeroIndex, i))); pause(0.1);
% %         hold off
% %     end
%     V(i) = fg.s;
%     w1s(i) = fg.w1;
%     w2s(i) = fg.w2;
% end

cftRegion = 200:500;
w1 = mean(w1s(cftRegion));
w2 = mean(w2s(cftRegion));
epsilon = 1e-3;
scatter(log(getScaledVariable(data.t(cftRegion), epsilon, L, 1)), V(cftRegion));
hl = lsline;
B = [(hl.YData(1) - hl.YData(2))/(hl.XData(1) - hl.XData(2)), (hl.YData(2)*hl.XData(1) - hl.YData(1)*hl.XData(2))/(hl.XData(1)-hl.XData(2))];
disp(strcat('1 = ', num2str(B(1) * (2 * pi)^2 * n)));
nuc = exp(4*pi^2 * 3 * B(2)/K);

load('rn3.mat');
rn = rn3;
[t, sortInds] = sort(t);
for i = 1:length(sortInds)
    rn(i, :) = rn3(sortInds(i), :);
end
for i = 1:length(t)
    r0(i) = sum(exp(-0i.*data.alphas).*abs(rn(i, :)/rn(i, zeroIndex))./(2* pi) * 1e-2);
    r1(i) = sum(exp(-1i.*data.alphas).*abs(rn(i, :)/rn(i, zeroIndex))./(2* pi) * 1e-2);
    r2(i) = sum(exp(-2i.*data.alphas).*abs(rn(i, :)/rn(i, zeroIndex))./(2* pi) * 1e-2);
end

plot(t, r1);
hold on
plot(t, chargeNegT(t, [epsilon nuc w1 w2], [l n K 1]));
