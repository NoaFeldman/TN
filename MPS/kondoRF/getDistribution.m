function distr = getDistribution(NRG, AV, dotVOp, dotCOp, leadOp, option)
    % For Wilson-local operator op, get its distribution along wilson
    % chain.
    if nargin == 5
        option = 'distribution';
    end
    
    if strcmp(option, 'distribution')
        sites = 2:length(NRG) - 1;
    else
        sites = length(NRG) - 1;
    end
    
    chainStart = contract(contract(AV, 3, dotVOp, 2), '12', AV, '12*');
    chainStart.info.itags = {'*', ''};
    cLink = contract(contract(NRG(1).AK, 3, dotCOp, 2), '3', NRG(1).AK, '3*');
    chainStart = contract(chainStart, '12', cLink, '13');
    for i = sites
        currChain = chainStart;
        for j = i:length(NRG) - 1
            if j == i
                %  |    __    |
                % |_|--|op|--|_|
                %  |          |
                %
                temp = contract(NRG(i).AK, 3, leadOp, 2);
                currChain = contract(contract(currChain, 1, temp, 1), '13', NRG(j).AK, '13*');
            else
                %  |      |
                % |_|----|_|
                %  |      |
                %
                currChain = contract(contract(currChain, 1, NRG(j).AK, 1), '13', NRG(j).AK, '13*');
            end
        end
        chainStart = contract(contract(chainStart, 1, NRG(i).AK, 1), '13', NRG(i).AK, '13*');
        currChain = contract(currChain, '12', getIdentity(currChain, 1), '12');
        if (isempty(currChain))
            distr(i) = 0;
        else
            distr(i) = getscalar(currChain);
        end
    end
end

            
            