function transformed = transform2(copies)
    transformed = QSpace(length(copies));
    [S,~]=getLocalSpace('Spin',0.5,'-A');
    NUp = contract(S(2)', '23', S(2), '13')*2;
    NDown = contract(S(3)', '23', S(3), '13')*2;
    for j = 1:length(copies)
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
        projUpUp = NUp;
        projUpUp.Q = [projUpUp.Q projUpUp.Q];
        projUpUp.info.itags = {strcat('A', int2str(j), 's'), ...
                               strcat('A', int2str(j), 's*'), ...
                               strcat('B', int2str(j), 's'), ...
                               strcat('B', int2str(j), 's*')};
        projDownDown = NDown;
        projDownDown.Q = [projDownDown.Q projDownDown.Q];
        projDownDown.info.itags = projUpUp.info.itags;
       upUp = contract(copies(j), '25', projUpUp, '24', [1 5 2 3 6 4]);
%        for i = 1:length(upUp.data)
%            upUp.daia{i} = upUp.daia{i} * (-1)^(upUp.Q{1}(i);
%        end
       transformed(j) = transformed(j) + upUp + ...
           contract(copies(j), '25', projDownDown, '24', [1 5 2 3 6 4]);
        
        % --O--                --O--                --O--
        %   |                    |                    |1
        %        >>  -1/sqrt(2)         + 1/sqrt(2)
        % --X--                --X--                --X--
        %   |                    |                    |0
        projDownUp = projUpUp;
        projDownUp.Q = [NDown.Q NUp.Q];
        downUp = contract(copies(j), '25', projDownUp, '24', [1 5 2 3 6 4]);
        upDown = contract(contract(downUp, 2, S(3), 2), '47', S(3)', '23', [1 5 2 3 6 4]);
        transformed(j) = transformed(j) + 1/sqrt(2)*(downUp - upDown);
        % --X--                --X--                --X--
        %   |                    |                    |0
        %        >>  1/sqrt(2)         + 1/sqrt(2)
        % --O--                --O--                --X--
        %   |                    |                    |1
        projUpDown = projUpUp;
        projUpDown.Q = [NUp.Q NDown.Q];
        upDown = contract(copies(j), '25', projUpDown, '24', [1 5 2 3 6 4]);
        downUp = contract(contract(upDown, 2, S(3)', 2), '47', S(3), '23', [1 5 2 3 6 4]);
        for i = 1:length(downUp.data)
           downUp.data{i} = downUp.data{i} * (-1)^(downUp.Q{1}(i));
        end
        transformed(j) = transformed(j) + 1/sqrt(2)*(downUp - upDown);
    end
end