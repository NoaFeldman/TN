function H = getH(N, h, JPM, JZ)
    % Create three arrays of length N,  one with single site operators, one
    % with operators meant to be paired with the left site and one with
    % operators meant to be paired with the right site. The right operators
    % always carry the JPM coefficient.
    
    [S,IS]=getLocalSpace('Spin',0.5,'-A');
    for i = 1 : length(S)
        S(i).info.itags = {'s*', 's', 'm*'};
    end
    
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
    
    H.single = QSpace(N);
    H.l2r = QSpace(N);
    H.r2l = QSpace(N);
    H.identity = QSpace(N);
    
    for i = 1 : N
        single = SZ_single;
        r2l = JZ * S(1) + JPM * S(2) + JPM * S(3);
        l2r = permute(conj(S(1)) + conj(S(2)) + conj(S(3)), [2 1 3]);
        id = identity;
        for j = 1:2
            single.info.itags(j) = strcat(int2str(i), single.info.itags(j));
            id.info.itags(j) = strcat(int2str(i), id.info.itags(j));
            r2l.info.itags(j) = strcat(int2str(i), r2l.info.itags(j));
            l2r.info.itags(j) = strcat(int2str(i), l2r.info.itags(j));
        end
        l2r.info.itags(3) = strcat(int2str(i), l2r.info.itags(3));
        r2l.info.itags(3) = strcat(int2str(i-1), r2l.info.itags(3));
        H.single(i) = single * -1 * h;
        H.identity(i) = id;
        if (i ~= 1) 
            H.r2l(i) = r2l;
        end
        if (i ~= N)
            H.l2r(i) = l2r;
        end
    end
    
    
        
