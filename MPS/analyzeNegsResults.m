% different params for negativity n
% For EE - in separate figs)


load('analyzeNegsTResults');
plot(t, r0, 'color', [0 0 0.3]);
hold on
plot(t, r1, 'color', [0 0 0.3]);
plot(t, r2, 'color', [0 0 0.3]);
plot(t, r3, 'color', [0 0 0.3]);

l = 1000;
n = 3;
K = 1;

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
plot(t, chargeNegT(t, [mean(epss1(tRangeInds)) mean(nucs21(tRangeInds)) mean(nucs22(tRangeInds)) mean(ws1(tRangeInds))], ...
    [l n K 0]), 'color', [0.9 0 0.4]);
plot(t, chargeNegT(t, [mean(epss1(tRangeInds)) mean(nucs21(tRangeInds)) mean(nucs22(tRangeInds)) mean(ws1(tRangeInds))], ...
    [l n K 1]), 'color', [0.9 0 0.4]);
plot(t, chargeNegT(t, [mean(epss1(tRangeInds)) mean(nucs21(tRangeInds)) mean(nucs22(tRangeInds)) mean(ws1(tRangeInds))], ...
    [l n K 2]), 'color', [0.9 0 0.4]);
plot(t, chargeNegT(t, [mean(epss1(tRangeInds)) mean(nucs21(tRangeInds)) mean(nucs22(tRangeInds)) mean(ws1(tRangeInds))], ...
    [l n K 3]), 'color', [0.9 0 0.4]);

hold off
plot(t, r1, 'color', [0 0 0.3]);
hold on
plot(t, chargeNegT(t, [epss1(:) nucs21(:) nucs22(:) ws1(:)], [l n K 1]));

hold off
