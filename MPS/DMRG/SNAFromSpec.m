function specResults = SNAFromSpec(dirName, firstStep, lastStep, stepT)
    stepsNum = lastStep - firstStep + 1;
    sz = -10:2:10;
    specResults.s = zeros(length(sz), stepsNum); 
    specResults.p = zeros(length(sz), stepsNum); 
    specResults.sFull = zeros(1, stepsNum); 
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
            dn = 1e-2;
            sNA = -(sum(val{i}.^(1+dn)) - sum(val{i}.^(1-dn))) / (2 * dn);
            if (~isempty(ind))
                specResults.s(ind, step+1 - firstStep) = sNA;
                specResults.p(ind, step+1 - firstStep) = pNA;
            end
            specResults.sFull(step+1 - firstStep) = specResults.sFull(step+1 - firstStep) + sNA;
            szForGaussianFit(i) = currSz/2;
            pForGaussianFit(i) = pNA;
            sForGaussianFit(i) = sNA;
        end
        if (length(val) >= 3)
            try
                f = fittype('sqrt(1/(2*pi*c))*exp(-(x-b)^2/(2*c))', 'independent', 'x', 'dependent', 'y');
                [fg, gof] = fit(szForGaussianFit.', real(pForGaussianFit).', f, 'StartPoint', [0.1 1]);
                specResults.sigmaN(step+1 - firstStep) = fg.c;
                specResults.avgN(step+1 - firstStep) = fg.b;
                specResults.gaussR(step+1 - firstStep) = gof.rsquare;
            catch exception
                specResults.avgN(step+1 - firstStep) = sum(szForGaussianFit.*pForGaussianFit);
            end
        end
    end
    specResults.t = firstStep:lastStep;
    specResults.t = specResults.t * stepT/2;
end