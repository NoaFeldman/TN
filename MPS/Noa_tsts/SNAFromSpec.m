function SNAFromSpec(dirName, firstStep, lastStep, stepT, figFileName, color)
    stepsNum = lastStep - firstStep + 1;
    sz = 0:2:10;
    s = zeros(length(sz), stepsNum); 
    p = zeros(length(sz), stepsNum); 
    sFull = zeros(1, stepsNum); 
    pFull = zeros(1, stepsNum); 
    for step = firstStep:lastStep
        file = load(strcat(dirName, '/step', int2str(step)));
        k = keys(file.spectrum);
        val = values(file.spectrum);
        for i = 1:length(val)
            currSz = str2num(k{i});
            ind = find(sz == currSz);
            pNA = sum(val{i});
            dn = 1e-2;
            sNA = -(sum(val{i}.^(1+dn)) - sum(val{i}.^(1-dn))) / (2 * dn);
            if (~isempty(ind))
                s(ind, step+1 - firstStep) = sNA;
                p(ind, step+1 - firstStep) = pNA;
            end
            sFull(step+1 - firstStep) = sFull(step+1 - firstStep) + sNA;
            pFull(step+1 - firstStep) = pFull(step+1 - firstStep) + pNA;
        end
    end
    t = firstStep:lastStep;
    t = t * stepT;
    plot(t, sFull, 'color', color);
    legendInfo{1} = ['sFull'];
    hold on
    for i = 1:length(sz)
        plot(t(:), s(i, :));
        legendInfo{i+1} = ['$2S^z = $' num2str(sz(i))];
%         plot(t(:), p(i, :));
%         legendInfo{i} = ['$2S^z = $' num2str(sz(i))];
    end
    set(gca, 'XScale', 'log');
%     legend(legendInfo, 'Interpreter', 'latex');
    xlabel('t', 'Interpreter', 'latex');
    ylabel('$S$', 'Interpreter', 'latex');
%     ylabel('$P(s^z)$', 'Interpreter', 'latex');
    savefig(figFileName);
    hold off
end