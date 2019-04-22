function distr = getDistribution(NRGPath, AV, AC, dotVOp, dotCOp, leadOp, T, option)
    % For Wilson-local operator op, get its distribution along wilson
    % chain.
    if nargin == 7
        option = 'distribution';
    end
    
        site = 0;
    NRG = struct('AK',{},'AT',{},'HK',{},'HT',{});
    while exist(strcat(NRGPath, '_', sprintf('%02d', site), '.mat')) == 2
        tmp = load(strcat(NRGPath, '_', sprintf('%02d', site)));
        NRG(site + 1).AK = tmp.AK; NRG(site + 1).AT = tmp.AT; NRG(site + 1).HK = tmp.HK; NRG(site + 1).HT = tmp.HT; 
        site = site + 1;
    end
    
    if strcmp(option, 'distribution')
        sites = 2:length(NRG) - 1;
    else
        sites = 2;
    end
    
    chainStart = contract(contract(AV, 3, dotVOp, 2), '12', AV, '12*');
    chainStart.info.itags = {'*', ''};
    cLink = contract(contract(AC, 3, dotCOp, 2), '3', AC, '3*');
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
        currChain = contract(contract(currChain, 1, NRG(end).AT, 1), '13', NRG(end).AT, '13*');
        rho0 = getThermalState(QSpace(NRG(end).HT), T);             
        distr(i) = getscalar(contract(currChain, '12', rho0, '12'));
        
        chainStart = contract(contract(chainStart, 1, NRG(i).AK, 1), '13', NRG(i).AK, '13*');
    end
end

            
            