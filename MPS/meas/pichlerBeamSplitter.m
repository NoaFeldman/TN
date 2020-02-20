function transformed = pichlerBeamSplitter(mergedCopies)
    % Eq. (2) in https://arxiv.org/pdf/1302.1187.pdf
    [S,IS]=getLocalSpace('Spin',0.5,'-A');
    NUp = contract(S(2)', '23', S(2), '13')*2;
    NDown = contract(S(3)', '23', S(3), '13')*2;
    splitter = QSpace();
    % --O--       --O--
    %   |           |
    %        >> 
    % --O--       --O--
    %   |           |
    
    % --X--       --X--
    %   |           |
    %        >> 
    % --X--       --X--
    %   |           |   
    splitter = splitter + ...
        contract(getIdentity(NUp, 1, NUp, 1), 3, getIdentity(NUp, 1, NUp, 1), '3*');
    splitter = splitter + ...
         contract(getIdentity(NDown, 1, NDown, 1), 3, getIdentity(NDown, 1, NDown, 1), '3*');
    % --O--                --O--                --O--
    %   |                    |                    |1
    %        >>  -1/sqrt(2)         + 1/sqrt(2)
    % --X--                --X--                --X--
    %   |                    |                    |0
    splitter = splitter - ...
        1/sqrt(2) * contract(getIdentity(NDown, 1, NUp, 1), 3, getIdentity(NDown, 1, NUp, 1), '3*');
    splitter = splitter + ...
        1/sqrt(2) * contract(getIdentity(NDown, 1, NUp, 1), 3, getIdentity(NUp, 1, NDown, 1), '3*');
    % --X--                --X--                --X--
    %   |                    |                    |0
    %        >>  1/sqrt(2)         + 1/sqrt(2)
    % --O--                --O--                --X--
    %   |                    |                    |1
    splitter = splitter + ...
        1/sqrt(2) * contract(getIdentity(NUp, 1, NDown, 1), 3, getIdentity( NUp, 1, NDown, 1), '3*');
    splitter = splitter + ...
        1/sqrt(2) * contract(getIdentity(NUp, 1, NDown, 1), 3, getIdentity(NDown, 1, NUp, 1), '3*');
    for i = 1:length(mergedCopies)
        transformed(i) = contract(mergedCopies(i), '23', splitter, '34', [1 3 4 2]);
        transformed(i).info.itags(2:3) = mergedCopies(i).info.itags(2:3);
    end
end