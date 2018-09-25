function [N1, N2] = getNegativityNs(psi, u1, v1, v2)
    % Get N matrices based on https://arxiv.org/pdf/1605.00674.pdf
    % psi is left canonical.
    u2 = v1 + 1;
    k = length(psi);
    while k > v2
        [psi, k] = shiftWorkingSite(psi, k, '<<');
    end
    E2 = getNegativityE(psi, u2, v2);
    N2 = getN(E2, 2);
    % N2 is now canonicalized such that it is the middle of the chain, move
    % working site to be v1.
    M = contract(psi(v1), 3, N2, 1);
    [psi(v1), N2] = orthoQS(M, [1, 2], '<<');
    N2.info.itags = {strcat(int2str(v1),'a*') N2.info.itags{2}, N2.info.itags{3}};
    N2 = QSpace(N2);
    psi(v1).info.itags = {psi(v1).info.itags{1}, psi(v1).info.itags{2}, strcat(int2str(v1),'a')};
    E1 = getNegativityE(psi, u1, v1);
    N1 = getN(E1, 1);
end

function N = getN(E, sysNum)
    % Contract two bottom legs, so now the charge on this leg is the charge
    % in E.
%     E = contract(E, '34', getIdentity(E, 3, E, 4), '12');
%     EData = E.data;
%     EQ = E.Q;
%     
%     NQ = {zeros(length(EData) ,1), zeros(length(EData) ,1), zeros(length(EData) ,1)};
%     NData = {};
%     
%     for i = 1:length(EData)
%         d1 = length(EData{i}(:, 1, 1));
%         d2 = length(EData{i}(1, :, 1));
%         d3 = length(EData{i}(1, 1, :));
%         [U, S, V] = svd(reshape(EData{i}, d1*d2, d3));
% %         if (d3 > d1*d2)
% %             S = S(:, 1:d1*d2);
% %         end
%         NMat = U*sqrt(S);
%         NMat = reshape(NMat, d1, d2, length(NMat(1, :)));
%         for j = 1:2
%             NQ{j}(i) = E.Q{j}(i);
%         end
%         % Flip sA direction
%         NQ{3}(i) = -E.Q{3}(i);
%         NData{i} = NMat;
%     end
%     N = QSpace;
%     N.info.itags = {E.info.itags{1}, E.info.itags{2}, strcat('sA', int2str(sysNum))};
%     N.info.qtype = '';
%     N.info.otype = '';
%     N.Q = NQ;
%     N.data = NData;
    [U, S, Vd, I] = svdQS(E, [3, 4], 'Nkeep', 256);
    I.svd2tr
    sqrtS = S;
    for i = 1:length(S.data)
        sqrtS.data{i} = zeros(length(S.data{i}));
        for j = 1:length(S.data{i})
            sqrtS.data{i}(j, j) = sqrt(S.data{i}(j));
        end
    end
    N = contract(QSpace(U), 3, QSpace(sqrtS), 1);
    N = contract(N, 3, getIdentity(N, 3, '-0', strcat('sA', int2str(sysNum))), 1);
%     if (sysNum == 2)
%         N = contract(N, 3, getIdentity(N, 3, '-0', strcat('sA', int2str(sysNum))), 1);
%     else
%         N.info.itags{3} = strcat('sA', int2str(sysNum), '*');
%     end
end