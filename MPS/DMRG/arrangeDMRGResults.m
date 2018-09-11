function specResults = arrangeDMRGResults(dirName, firstStep, lastStep, stepT)
    stepsNum = lastStep - firstStep + 1;
    sz = -10:2:10;
    specResults.alphas = -3.14:0.01:3.14;
    specResults.s = zeros(length(sz), stepsNum); 
    specResults.s1Charge = zeros(length(sz), stepsNum); 
    specResults.sFull = zeros(1, stepsNum);
    specResults.s1Alpha = zeros(length(specResults.alphas), stepsNum); 
    specResults.sigmaN = zeros(1, stepsNum);
    specResults.avgN = zeros(1, stepsNum);
    specResults.gaussR = zeros(1, stepsNum);
    for step = firstStep:lastStep
        file = load(strcat(dirName, '/step', int2str(step)));
        k = keys(file.spectrum);
        val = values(file.spectrum);
        szForGaussianFit = zeros(1, length(val)); 
        pForGaussianFit = zeros(1, length(val));
        sForGaussianFit = zeros(1, length(val));
        for i = 1:length(val)
            currSz = str2num(k{i});
            ind = find(sz == currSz);
            pNA = sum(val{i});
%             dn = 1e-2;
%             sNA = -(sum(val{i}.^(1+dn)) - sum(val{i}.^(1-dn))) / (2 * dn);
            sNA = -sum(val{i}.*log(val{i}));
            if (~isempty(ind))
                specResults.s(ind, step+1 - firstStep) = sNA;
                specResults.s1Charge(ind, step+1 - firstStep) = pNA;
            end
            specResults.sFull(step+1 - firstStep) = specResults.sFull(step+1 - firstStep) + sNA;
            szForGaussianFit(i) = currSz/2;
            pForGaussianFit(i) = pNA;
            sForGaussianFit(i) = sNA;
        end
        specResults.s1Alpha(:, step+1 - firstStep) = discreteFourier(specResults.alphas, sz/2, specResults.s1Charge(:, step+1 - firstStep).');
        if (length(val) >= 3)
            try
                f = fittype('sqrt(1/(2*pi*c))*exp(-(x)^2/(2*c))', 'independent', 'x', 'dependent', 'y');
                [fg, gof] = fit(szForGaussianFit.', real(pForGaussianFit).', f, 'StartPoint', [1]);
                specResults.sigmaN(step+1 - firstStep) = fg.c;
%                 specResults.avgN(step+1 - firstStep) = fg.b;
                specResults.gaussR(step+1 - firstStep) = gof.rsquare;
                plot(fg, szForGaussianFit.', real(pForGaussianFit).'); pause(0.1);
            catch exception
%                 specResults.avgN(step+1 - firstStep) = sum(szForGauspecssianFit.*pForGaussianFit);
            end
        end
    end
    specResults.t = firstStep:lastStep;
    specResults.t = specResults.t * stepT/2;
end

function sAlpha = discreteFourier(alphas, sz, sCharge)
    sAlpha = sum(exp(1i.*alphas.*sz.').*sCharge.', 1);
end