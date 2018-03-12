function SNAFromSpec(dirName, stepsNum, stepT, figFileName)
    s0 = zeros(1, stepsNum+1); 
    for step = 0:stepsNum
        s = load(strcat(dirName, '/step', int2str(step)));
        val = values(s.spectrum);
        for i = 1:length(val)
            s0(step+1) = s0(step+1) - sum(val{i}.*log2(val{i}));
    %         dn = 1e-2;
    %         s0(step+1) = -(sum(spec0{i}.^(1+dn)) - sum(spec0{i}.^(1-dn))) / (2 * dn);
        end
    end
    t = 0:stepsNum;
    t = t * stepT;
    scatter(t(:), s0(:));
    xlabel('t', 'Interpreter', 'latex');
    ylabel('$S(0)$', 'Interpreter', 'latex');
    savefig(figFileName);
    hold off
end