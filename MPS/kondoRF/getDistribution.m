function dist = getDistribution(NRGD, AV, dotVOp, dotCOp, leadOp)
    % For Wilson-local operator op, get its distribution along wilson
    % chain.
    
    % first site is the QD conductance level - we skip it.
    chainStart = contract(contract(AV, 1, dotVOp, 2), 2, AV, '1*');
    cLink = contract(contract(NRGD(1).AK, 3, dotCOp, 2), '3', NRGD(1).AK, '3*');
    chainStart = contract(chainStart, '12', cLink, '13');
    for i = 2:length(NRGD) - 1
        currChain = chainStart;
        for j = 2:length(NRGD) - 1
            if j == i
                %  |    __    |
                % |_|--|op|--|_|
                %  |          |
                %
                newLink = contract(contract(NRGD(i).AK, 3, leadOp, 2), '3', NRGD(i).AK, '3*');
            else
                %  |      |
                % |_|----|_|
                %  |      |
                %
                newLink = contract(NRGD(j).AK, '3', NRGD(j).AK, '3*');
            end
            currChain = contract(currChain, '12', newLink, '13');
        end
        currChain = contract(currChain, '12', getIdentity(currChain, 1), '12');
        if (isempty(currChain))
            dist(i) = 0;
        else
            dist(i) = currChain.data{1};
        end
    end
end

            
            