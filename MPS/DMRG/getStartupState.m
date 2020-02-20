function psi = getStartupState(N, d, m, bc)
    % Create a product state (up-down-up-down-up...) for n sites with spin m.
    % Assuming n even.
    % Assuming m >= 0.
    % Bring the state to a left canonical form (except for the rightmost site).
    if nargin == 1
        d = 2
        m = 0;
        bc = 'open';
    end
    if nargin == 2
        m = 0;
        bc = 'open';
    end
    if d == 2
        if strcmp(bc, 'open')
            fullState1 = QSpace();
            fullState1.info.qtype = '';
            fullState1.info.otype = '';
            fullState1.Q{1} = 0;
            fullState1.info.itags{1} = '0a';
            for i = 2:N+1
                if (mod(i, 2) == 0)
                    fullState1.Q{i} = 1;
                else
                    fullState1.Q{i} = -1;
                end
                fullState1.info.itags{i} = strcat(int2str(i-1), 's');
            end
            fullState1.Q{N+2} = 0;
            fullState1.info.itags{N+2} = strcat(int2str(N), 'a');
            fullState1.data{1} = 1;

            fullStateFlipped1 = fullState1;
            for i = 2:N+1
                fullStateFlipped1.Q{i} = -fullState1.Q{i};
            end
            fullState1 = fullState1 + fullStateFlipped1;

            fullState2 = QSpace();
            fullState2.info.qtype = '';
            fullState2.info.otype = '';
            fullState2.Q{1} = 0;
            fullState2.info.itags{1} = '0a';
            for i = 2:N+1
                if (mod(i, 4) == 0 || mod(i+1, 4) == 0)
                    fullState2.Q{i} = 1;
                else
                    fullState2.Q{i} = -1;
                end
                fullState2.info.itags{i} = strcat(int2str(i-1), 's');
            end
            fullState2.Q{N+2} = 0;
            fullState2.info.itags{N+2} = strcat(int2str(N), 'a');
            fullState2.data{1} = 1;

            fullStateFlipped2 = fullState2;
            for i = 2:N+1
                fullStateFlipped2.Q{i} = -fullState2.Q{i};
            end
            fullState2 = fullState2 + fullStateFlipped2;

            fullState3 = QSpace();
            fullState3.info.qtype = '';
            fullState3.info.otype = '';
            fullState3.Q{1} = 0;
            fullState3.info.itags{1} = '0a';
            for i = 2:N+1
                if (mod(i-1, 4) == 0 || mod(i, 4) == 0)
                    fullState3.Q{i} = 1;
                else
                    fullState3.Q{i} = -1;
                end
                fullState3.info.itags{i} = strcat(int2str(i-1), 's');
            end
            fullState3.Q{N+2} = 0;
            fullState3.info.itags{N+2} = strcat(int2str(N), 'a');
            fullState3.data{1} = 1;

            fullStateFlipped3 = fullState3;
            for i = 2:N+1
                fullStateFlipped3.Q{i} = -fullState3.Q{i};
            end
            fullState3 = fullState3 + fullStateFlipped3;

    %         fullState = fullState1 + fullState2 + fullState3;
            fullState = fullState3 - fullState1;
            psi = QSpace(N);
            for i = 1:N-1
                [psi(i), fullState, ~] = orthoQS(fullState, [1 2], '>>', ...
                                'Nkeep', 1024, 'itag', strcat(int2str(i), 'a'));
            end
            psi(N) = fullState;
            norm = getOverlap(psi, psi);
            for i = 1:length(psi(N).data)
                psi(N).data{i} = psi(N).data{i}/sqrt(norm);
            end

            if m == -1
                % for now only working for |m|, if we go to the asymmetric case
                % fix.
                [S,IS]=getLocalSpace('Spin',0.5,'-A');
                S(2).Q{3} = -S(2).Q{3};
                S(2).info.itags{3} = '';
                S(2).info.itags = {strcat(int2str(N), 's'), strcat(int2str(N), 's*'), 'temp'};
                temp = contract(psi(N), 2, S(2), 2);
                id = getIdentity(temp, 2, temp, 4, psi(N).info.itags{3});
                psi(N) = contract(temp, '24', id, '12*')
            end
            if m ~= 0 && m ~= -1
                disp('------------------------Unsupported m! Can be fixed in getStartupState.------------------');
            end
        else

            fullState1 = QSpace();
            fullState1.info.qtype = '';
            fullState1.info.otype = '';
            fullState1.Q{1} = 0;
            fullState1.info.itags{1} = '0a';
            for i = 1:N
                fullState1.Q{2*i} = getQ(i, N, 'u');
                fullState1.Q{2*i+1} = getQ(i, N, 'd');
                fullState1.info.itags{2*i} = strcat(int2str(i), 'sUp');
                fullState1.info.itags{2*i+1} = strcat(int2str(i), 'sDown');
            end
            fullState1.Q{2*N+2} = 0;
            fullState1.info.itags{2*N+2} = strcat(int2str(N), 'a');
            fullState1.data{1} = 1;

            fullStateFlipped1 = fullState1;
            for i = 2:2*N+1
                fullStateFlipped1.Q{i} = -fullState1.Q{i};
            end
    %         fullState1 = fullState1 + fullStateFlipped1;

            fullState2 = fullState1;
            up1ind = find(strcmp(fullState1.info.itags, '1sUp'));
            up3ind = find(strcmp(fullState1.info.itags, '3sUp'));
            fullState2.Q{up1ind} = -fullState2.Q{up1ind};
            fullState2.Q{up3ind} = -fullState2.Q{up3ind};

            fullState3 = QSpace();
            fullState3.info.qtype = '';
            fullState3.info.otype = '';
            fullState3.Q{1} = 0;
            fullState3.info.itags{1} = '0a';
            for i = 1:N
                if (mod(i, 2) == 0)
                    fullState3.Q{2*i} = 1;
                    fullState3.Q{2*i+1} = -1;
                else
                    fullState3.Q{2*i} = -1;
                    fullState3.Q{2*i+1} = 1;
                end
                fullState3.info.itags{2*i} = strcat(int2str(i), 'sUp');
                fullState3.info.itags{2*i+1} = strcat(int2str(i), 'sDown');
            end
            fullState3.Q{2*N+2} = 0;
            fullState3.info.itags{2*N+2} = strcat(int2str(N), 'a');
            fullState3.data{1} = 1;

            fullState4 = fullState3;
            up1ind = find(strcmp(fullState1.info.itags, '1sUp'));
            up2ind = find(strcmp(fullState1.info.itags, '2sUp'));
            fullState4.Q{up1ind} = -fullState4.Q{up1ind};
            fullState4.Q{up2ind} = -fullState4.Q{up2ind};

    %         fullStateFlipped3 = fullState3;
    %         for i = 2:2*N+1
    %             fullStateFlipped3.Q{i} = -fullState3.Q{i};
    %         end
    %         fullState3 = fullState3 + fullStateFlipped3;

            fullState = fullState1 + fullState2 + fullState3 + fullState4;
            psi = QSpace(N);
            for i = 1:N-1
                [psi(i), fullState, ~] = orthoQS(fullState, [1 2 3], '>>', ...
                                'Nkeep', 1024, 'itag', strcat(int2str(i), 'a'));
            end
            psi(N) = fullState;
            norm = getOverlap(psi, psi);
            for i = 1:length(psi(N).data)
                psi(N).data{i} = psi(N).data{i}/sqrt(norm);
            end

            if m == -1
                % for now only working for |m|, if we go to the asymmetric case
                % fix.
                [S,IS]=getLocalSpace('Spin',0.5,'-A');
                S(2).Q{3} = -S(2).Q{3};
                S(2).info.itags{3} = '';
                S(2).info.itags = {strcat(int2str(N), 'sDown'), strcat(int2str(N), 'sDown*'), 'temp'};
                temp = contract(psi(N), 3, S(2), 2);
                id = getIdentity(temp, 3, temp, 5, psi(N).info.itags{3});
                psi(N) = contract(temp, '35', id, '12*');
            end
            if m ~= 0 && m ~= -1
                disp('------------------------Unsupported m! Can be fixed in getStartupState.------------------');
            end
        end
    end
    if d == 3
        fullState1 = QSpace();
        fullState1.info.qtype = '';
        fullState1.info.otype = '';
        fullState1.Q{1} = 0;
        fullState1.info.itags{1} = '0a';
        for i = 2:N+1
            if (mod(i, 4) == 0)
                fullState1.Q{i} = 2;
            elseif ((mod(i, 4) == 1) || (mod(i, 4) == 3))
                fullState1.Q{i} = 0;
            else
                fullState1.Q{i} = -2;
            end
            fullState1.info.itags{i} = strcat(int2str(i-1), 's');
        end
        fullState1.Q{N+2} = 0;
        fullState1.info.itags{N+2} = strcat(int2str(N), 'a');
        fullState1.data{1} = 1;
        
        fullState2 = QSpace();
        fullState2.info.qtype = '';
        fullState2.info.otype = '';
        fullState2.Q{1} = 0;
        fullState2.info.itags{1} = '0a';
        for i = 2:N+1
            if (mod(i, 4) == 1)
                fullState2.Q{i} = 2;
            elseif ((mod(i, 4) == 0) || (mod(i, 4) == 2))
                fullState2.Q{i} = 0;
            else
                fullState2.Q{i} = -2;
            end
            fullState2.info.itags{i} = strcat(int2str(i-1), 's');
        end
        fullState2.Q{N+2} = 0;
        fullState2.info.itags{N+2} = strcat(int2str(N), 'a');
        fullState2.data{1} = 1;
        
        fullState = fullState2 - fullState1;
        psi = QSpace(N);
        for i = 1:N-1
            [psi(i), fullState, ~] = orthoQS(fullState, [1 2], '>>', ...
                            'Nkeep', 1024, 'itag', strcat(int2str(i), 'a'));
        end
        psi(N) = fullState;
        norm = getOverlap(psi, psi);
        for i = 1:length(psi(N).data)
            psi(N).data{i} = psi(N).data{i}/sqrt(norm);
        end
    end
end

function q = getQ(ind, N, chainHalf)
    if strcmp(chainHalf, 'u')
        ind = N*2 - ind + 1;
    end
    if (mod(ind, 4) == 0 || mod(ind+1, 4) == 0)
        q = 1;
    else
        q = -1;
    end
end