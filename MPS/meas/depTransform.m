function transformed = depTransform(psi, copyNum)
% not the same for both copies - hve to go for expression on paper
% merged should be particle number conserving
    omega = exp(2 * pi * 1i / copyNum);
    [S,IS]=getLocalSpace('Spin',0.5,'-A');
    projUp = contract(S(2)', '23', S(2), '13')*2;
    aDagg = S(3);
    transformed = QSpace(length(psi));

    for i = 1:length(psi)
        Ni = particleNumber(psi, i, projUp);
        if i == 1
            transformed(i) = emptySite(0, i, copyNum);
        else
            transformed(i) = emptySite(transformed(i-1).Q{copyNum + 2}, i, copyNum);
            if i == length(psi)
                transformed(i).info.itags{copyNum + 2} = strcat(int2str(i), 'a');
                transformed(i).Q{copyNum + 2} = -transformed(i).Q{copyNum + 2};
            end
        end
        for n = 1:copyNum
            m1Term = QSpace();
            for k = 1:copyNum
                curr = contract(transformed(i), k+1, aDagg, 2, [1:k, copyNum + 2, k+1:copyNum + 1, copyNum + 3]);
                if i == length(psi)
                    curr.info.itags{copyNum + 3} =  '';
                    curr.Q{copyNum + 3} = -curr.Q{copyNum + 3};
                end
                curr.info.itags{k+1} = transformed(i).info.itags{k+1};
                indStr = strcat(int2str(copyNum+2), int2str(copyNum+3));
                if i == length(psi)
                    curr = contract(curr, indStr, ...
                        getIdentity(real(curr), copyNum+2, real(curr), copyNum+3, curr.info.itags{copyNum+2}), '12*');
                else
                    curr = contract(curr, indStr, ...
                        getIdentity(real(curr), copyNum+2, real(curr), copyNum+3, curr.info.itags{copyNum+2}), '12');
                end
                m1Term = m1Term + omega^(n*k)*curr;
            end
            transformed(i) = m1Term * Ni + transformed(i)*(1 - Ni);
        end
    end
end

function merged = merge(copy1, copy2)
    idIn = getIdentity(real(copy1), 1, real(copy2), 1, copy1.info.itags{1});
    idOut = getIdentity(real(copy1), 3, real(copy2), 3, copy1.info.itags{3});
    merged = contract(contract(idIn, '1*', copy1, 1), 1, copy2, 1);
    if ~isempty(strfind(merged.info.itags{3}, '*')) % Mid chain (for OC on the right)
        merged = contract(merged, '35', ...
            idOut, '12');
    else
        merged = contract(merged, '35', ...
            idOut, '12*');
    end
    merged.info.itags{2} = strcat('1', merged.info.itags{2});
    merged.info.itags{3} = strcat('2', merged.info.itags{3});
end

function qspace = emptySite(sIn, siteInd, copyNum)
    qspace = QSpace();
    qspace.Q = {sIn};
    for j = 1:copyNum
        qspace.Q{j + 1} = -1.*ones(length(sIn), 1);
    end
    qspace.Q{copyNum + 2} = sIn - copyNum;
    for j = 1:length(sIn)
        qspace.data{j} = 1/sqrt(length(sIn)); % Or maybe 1?????????????????????????????
    end
    qspace.info.itags = {strcat(int2str(siteInd - 1), 'a')};
    for j = 1:copyNum
        qspace.info.itags{j + 1} = strcat(int2str(j), '_', int2str(siteInd), 's');
    end
	qspace.info.itags{copyNum + 2} = strcat(int2str(siteInd), 'a*');
    qspace.info.qtype = '';
    qspace.info.otype = '';   
end

function state = getVacState(N, S)
    state = QSpace(N);
    m = 0;
    for i = 1:N
        state(i).Q = {m -1 (m-1)};
        m = m - 1;
        state(i).data = {1};
        tag1 = strcat(int2str(i - 1), 'a');
        tag2 = strcat(int2str(i), 's');
        tag3 = strcat(int2str(i), 'a', '*');
        state(i).info.itags = {tag1, tag2, tag3};
        state(i).info.qtype = '';
        state(i).info.otype = '';
    end
end

function Ni = particleNumber(psi, i, NOp)
    psi(i) = contract(psi(i), 2, NOp, 2, [1 3 2]);
    curr = contract(psi(1), '12', psi(1), '12*', [2 1]);
    for j = 2:length(psi)-1
        curr = contract(contract(curr, 2, psi(j), 1), '12', psi(j), '12*', [2 1]);
    end
    curr = contract(contract(curr, 2, psi(length(psi)), 1), '123', psi(length(psi)), '123*');
    Ni = getscalar(curr);
end