% different params for negativity n
% For EE - in separate figs)
load('analyzeNegsTResults');

l = 1000;
n = 3;
K = 1;

%todo eps = 30!
tic;
for i = 1:10:length(t)
    for j = 1:10:length(t)
        for k = 1:10:length(t)
            plot(t, chargeNegT(t, [30 nucs1(i) nucs1(j) ws1(k)], [l n K 1]));
            hold on
            toc;
        end
    end
end
plot(t, r1, 'color', 'c');
hold off

plot(t, r0, 'color', [0 0 0.3]);
hold on
plot(t, r1, 'color', [0 0 0.3]);
plot(t, r2, 'color', [0 0 0.3]);
plot(t, r3, 'color', [0 0 0.3]);

tRangeInds = 10:20;
% For one scaling constant
plot(t, chargeNegT(t, [mean(epss1(tRangeInds)) mean(nucs1(tRangeInds)) mean(nucs1(tRangeInds)) mean(ws1(tRangeInds))], ...
    [l n K 0]), 'color', [0 0.9 0.4]);
plot(t, chargeNegT(t, [mean(epss1(tRangeInds)) mean(nucs1(tRangeInds)) mean(nucs1(tRangeInds)) mean(ws1(tRangeInds))], ...
    [l n K 1]), 'color', [0 0.9 0.4]);
plot(t, chargeNegT(t, [mean(epss1(tRangeInds)) mean(nucs1(tRangeInds)) mean(nucs1(tRangeInds)) mean(ws1(tRangeInds))], ...
    [l n K 2]), 'color', [0 0.9 0.4]);
plot(t, chargeNegT(t, [mean(epss1(tRangeInds)) mean(nucs1(tRangeInds)) mean(nucs1(tRangeInds)) mean(ws1(tRangeInds))], ...
    [l n K 3]), 'color', [0 0.9 0.4]);

% Two scaling constants
plot(t, chargeNegT(t, [mean(epss2(tRangeInds)) mean(nucs21(tRangeInds)) mean(nucs22(tRangeInds)) mean(ws2(tRangeInds))], ...
    [l n K 0]), 'color', [0.9 0 0.4]);
plot(t, chargeNegT(t, [mean(epss2(tRangeInds)) mean(nucs21(tRangeInds)) mean(nucs22(tRangeInds)) mean(ws2(tRangeInds))], ...
    [l n K 1]), 'color', [0.9 0 0.4]);
plot(t, chargeNegT(t, [mean(epss2(tRangeInds)) mean(nucs21(tRangeInds)) mean(nucs22(tRangeInds)) mean(ws2(tRangeInds))], ...
    [l n K 2]), 'color', [0.9 0 0.4]);
plot(t, chargeNegT(t, [mean(epss2(tRangeInds)) mean(nucs21(tRangeInds)) mean(nucs22(tRangeInds)) mean(ws2(tRangeInds))], ...
    [l n K 3]), 'color', [0.9 0 0.4]);

hold off
