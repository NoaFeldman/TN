function H = getH(N, d, params, bc)
    % Create three arrays of length N,  one with single site operators, one
    % with operators meant to be paired with the left site and one with
    % operators meant to be paired with the right site. The right operators
    % always carry the JPM coefficient.
    if d == 2
        h = params(1);
        JPM = params(2);
        JZ = params(3);
        J2PM = params(4);
        J2Z = params(5);
        if strcmp(bc, 'open')
            H.single = QSpace(N);
            H.l2r = QSpace(N);
            H.r2l = QSpace(N);
            H.l2r2 = QSpace(N);
            H.r2l2 = QSpace(N);
            H.identity = QSpace(N);
            [S,IS]=getLocalSpace('Spin',0.5,'-A');
            for i = 1 : length(S)
                S(i).info.itags = {'s', 's*', 'm*'};
            end
    %         S(3).data{1} = -S(3).data{1}; % ?????????????????????????????????????????

            SZ_single  = S(1);
            % remove redundant 'm' leg
            SZ_single.Q = SZ_single.Q(1:2);
            SZ_single.info.itags = SZ_single.info.itags(1:2);
            % prepare the single site identity operator.
            identity = QSpace;
            identity.Q = {[1; -1], [1; -1]};
            identity.data = {1; 1};
            identity.info.qtype = '';
            identity.info.itags = SZ_single.info.itags;

            for i = 1 : N
                single = SZ_single;
                r2l = JZ * S(1) + JPM * S(2) + JPM * S(3);
                l2r = permute(conj(S(1)) + conj(S(2)) + conj(S(3)), [2 1 3]);
                r2l2 = J2Z * S(1) + J2PM * S(2) + J2PM * S(3);
                l2r2 = permute(conj(S(1)) + conj(S(2)) + conj(S(3)), [2 1 3]);
                id = identity;
                for j = 1:2
                    single.info.itags(j) = strcat(int2str(i), single.info.itags(j));
                    id.info.itags(j) = strcat(int2str(i), id.info.itags(j));
                    r2l.info.itags(j) = strcat(int2str(i), r2l.info.itags(j));
                    l2r.info.itags(j) = strcat(int2str(i), l2r.info.itags(j));
                    if ~isempty(r2l2)
                        r2l2.info.itags(j) = strcat(int2str(i), r2l2.info.itags(j));
                    end
                    l2r2.info.itags(j) = strcat(int2str(i), l2r2.info.itags(j));
                end
                l2r.info.itags(3) = strcat(int2str(i), l2r.info.itags(3));
                r2l.info.itags(3) = strcat(int2str(i-1), r2l.info.itags(3));
                if ~isempty(r2l2)
                    r2l2.info.itags(3) = strcat(int2str(i-2), r2l2.info.itags(3));
                end
                l2r2.info.itags(3) = strcat(int2str(i), l2r2.info.itags(3));
                H.single(i) = single * -1 * h;
                H.identity(i) = id;
                if (i ~= 1) 
                    H.r2l(i) = r2l;
                    if (i ~= 2)
                        H.r2l2(i) = r2l2;
                    end
                end
                if (i ~= N)
                    H.l2r(i) = l2r;
                    if (i ~= N-1)
                        H.l2r2(i) = l2r2;
                    end
                end
            end
        else
            H.single = QSpace(N);
            H.l2rUp = QSpace(N);
            H.l2rDown = QSpace(N);
            H.r2lUp = QSpace(N);
            H.r2lDown = QSpace(N);
            H.l2r2Up = QSpace(N);
            H.l2r2Down = QSpace(N);
            H.r2l2Up = QSpace(N);
            H.r2l2Down = QSpace(N);

            [S,IS]=getLocalSpace('Spin',0.5,'-A');
            % Add extra leg for combining up and down operators
            IS.E.Q = [IS.E.Q(1:2) {[0; 0]}];
            IS.E.info.itags = {'', '*', ''};
            S(2).Q = [S(2).Q(1:3) {0}];
            S(2).info.itags = {'', '*', '*', '*'};
            S(3).Q = [S(3).Q(1:3) {0}];
            S(3).info.itags = {'', '*', '*', '*'};
            S(3).data{1} = -S(3).data{1}; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            % operators for mid-chain
            single = (contract(S(1), 3, IS.E, 3) + contract(IS.E, 3, S(1), 3)) * -1 * h;
            for i = 1:N
                H.single(i) = single;
                H.single(i).info.itags = {strcat(int2str(i), 'sUp'), strcat(int2str(i), 'sUp*'), ...
                                        strcat(int2str(i), 'sDown'), strcat(int2str(i), 'sDown*')};
            end

            S(1).Q = [S(1).Q(1:3) {[0; 0]}];
            S(1).info.itags = {'', '*', '*', '*'};
            l2rUp = contract(JZ * S(1) + JPM * S(2) + JPM * S(3), 4, IS.E, 3, [1 2 4 5 3]);
            l2rDown = contract(IS.E, 3, JZ * S(1) + JPM * S(2) + JPM * S(3), 4);
            for i = 1:N-1
                H.l2rUp(i) = l2rUp;
                H.l2rUp(i).info.itags = {strcat(int2str(i), 'sUp'), strcat(int2str(i), 'sUp*'), ...
                                        strcat(int2str(i), 'sDown'), strcat(int2str(i), 'sDown*'), ...
                                        strcat(int2str(i), 'm*')};
                H.l2rDown(i) = l2rDown;
                H.l2rDown(i).info.itags = {strcat(int2str(i), 'sUp'), strcat(int2str(i), 'sUp*'), ...
                                        strcat(int2str(i), 'sDown'), strcat(int2str(i), 'sDown*'), ...
                                        strcat(int2str(i), 'm*')};
            end
            r2lUp = contract(S(1) + S(2) + S(3), 4, IS.E, 3, [1 2 4 5 3]);
            r2lUp.info.itags{5} = '';
            r2lUp.Q{5} = -r2lUp.Q{5};
            r2lDown = contract(IS.E, 3, S(1) + S(2) + S(3), 4);
            r2lDown.info.itags{5} = '';
            r2lDown.Q{5} = -r2lDown.Q{5};
            for i = 2:N
                H.r2lUp(i) = r2lUp;
                H.r2lUp(i).info.itags = {strcat(int2str(i), 'sUp'), strcat(int2str(i), 'sUp*'), ...
                                        strcat(int2str(i), 'sDown'), strcat(int2str(i), 'sDown*'), ...
                                        strcat(int2str(i-1), 'm')};
                H.r2lDown(i) = r2lDown;
                H.r2lDown(i).info.itags = {strcat(int2str(i), 'sUp'), strcat(int2str(i), 'sUp*'), ...
                                        strcat(int2str(i), 'sDown'), strcat(int2str(i), 'sDown*'), ...
                                        strcat(int2str(i-1), 'm')};
            end

            l2r2Up = contract(J2Z * S(1) + J2PM * S(2) + J2PM * S(3), 4, IS.E, 3, [1 2 4 5 3]);
            l2r2Down = contract(IS.E, 3, J2Z * S(1) + J2PM * S(2) + J2PM * S(3), 4);
            for i = 1:N-1
                H.l2r2Up(i) = l2r2Up;
                H.l2r2Up(i).info.itags = {strcat(int2str(i), 'sUp'), strcat(int2str(i), 'sUp*'), ...
                                        strcat(int2str(i), 'sDown'), strcat(int2str(i), 'sDown*'), ...
                                        strcat(int2str(i), 'm*')};
                H.l2r2Down(i) = l2r2Down;
                H.l2r2Down(i).info.itags = {strcat(int2str(i), 'sUp'), strcat(int2str(i), 'sUp*'), ...
                                        strcat(int2str(i), 'sDown'), strcat(int2str(i), 'sDown*'), ...
                                        strcat(int2str(i), 'm*')};
            end
            r2l2Up = contract(S(1) + S(2) + S(3), 4, IS.E, 3, [1 2 4 5 3]);
            r2l2Up.info.itags{5} = '';
            r2l2Up.Q{5} = -r2l2Up.Q{5};
            r2l2Down = contract(IS.E, 3, S(1) + S(2) + S(3), 4);
            r2l2Down.info.itags{5} = '';
            r2l2Down.Q{5} = -r2l2Down.Q{5};
            for i = 2:N
                H.r2l2Up(i) = r2l2Up;
                H.r2l2Up(i).info.itags = {strcat(int2str(i), 'sUp'), strcat(int2str(i), 'sUp*'), ...
                                        strcat(int2str(i), 'sDown'), strcat(int2str(i), 'sDown*'), ...
                                        strcat(int2str(i-2), 'm')};
                H.r2l2Down(i) = r2l2Down;
                H.r2l2Down(i).info.itags = {strcat(int2str(i), 'sUp'), strcat(int2str(i), 'sUp*'), ...
                                        strcat(int2str(i), 'sDown'), strcat(int2str(i), 'sDown*'), ...
                                        strcat(int2str(i-2), 'm')};
            end

            % Couple operators on edges (fold into each other)
            [S3,~]=getLocalSpace('Spin',0.5,'-A');
            S1Dagger = S3(1);
            S1Dagger.info.itags{3} = '';
            single = JPM * (contract(S3(2), 3, S3(2)', 3) + contract(S3(2)', 3, S3(2), 3)) + ...
                     JZ * contract(S3(1), 3, S1Dagger, 3);
            H.single(1) = H.single(1) + single;
            H.single(N) = H.single(N) + single;
            H.l2rUp(1) = H.l2rUp(1) + l2r2Down;
            H.l2rUp(N-1) = H.l2rUp(N-1) + l2r2Down;
            H.l2rDown(1) = H.l2rDown(1) + l2r2Up;
            H.l2rDown(N-1) = H.l2rDown(N-1) + l2r2Up;
        end
    end
    % AKLT
    if d == 3
        D = params(1);
        H.single = QSpace(N);
        H.l2r = QSpace(N);
        H.r2l = QSpace(N);
        
        H.l2r2 = QSpace(N);
        H.r2l2 = QSpace(N);
        
        [S, IS]=getLocalSpace('Spin', 1, '-A');
        for i = 1 : length(S)
            S(i).info.itags = {'s', 's*', 'm*'};
        end
        SZ_single  = S(1);
        % remove redundant 'm' leg
        SZ_single.Q = SZ_single.Q(1:2);
        SZ_single.info.itags = SZ_single.info.itags(1:2);
        % prepare the single site identity operator.

        for i = 1 : N
            single = SZ_single * SZ_single;
            r2l = S(1) + S(2) + S(3);
            l2r = permute(conj(S(1)) + conj(S(2)) + conj(S(3)), [2 1 3]);
            for j = 1:2
                single.info.itags(j) = strcat(int2str(i), single.info.itags(j));
                r2l.info.itags(j) = strcat(int2str(i), r2l.info.itags(j));
                l2r.info.itags(j) = strcat(int2str(i), l2r.info.itags(j));
            end
            l2r.info.itags(3) = strcat(int2str(i), l2r.info.itags(3));
            r2l.info.itags(3) = strcat(int2str(i-1), r2l.info.itags(3));
            H.single(i) = single * D;
            if (i ~= 1) 
                H.r2l(i) = r2l;
            end
            if (i ~= N)
                H.l2r(i) = l2r;
            end
        end
    end
end
    
        
