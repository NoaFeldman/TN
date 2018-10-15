function r = fluxRenyi3(cicj, alpha, u1, v1, v2)
    % Based on B6 in https://arxiv.org/pdf/1804.00632.pdf
    % Assuming A1 and A2are adjacent.
    u2 = v1+1;
    Gplus = getG(1, cicj, u1, v1, v2);
    Gminus = getG(-1, cicj, u1, v1, v2);
    N = eye(v2 - u1 + 1); 
    I = eye (v2 - u1 + 1);
    Gpp = (I + Gplus)/2;
    Gpp2 = Gpp^2;
    Gpp3 = Gpp2 * Gpp;
    Gpm = (I - Gplus)/2;
    Gpm2 = Gpm^2;
    Gpm3 = Gpm2 * Gpm;
    Tr = -1/2 * det(Gpm3 + expm(1i .* alpha .* N) * Gpp3) + ...
        3/2 * det(Gpm2 * (I - Gminus)/2 + expm(1i .* alpha .* N) * Gpp2 * (I + Gminus) / 2);
    
%     Tr = 0;
%     Gs = {Gminus, Gplus};
%     signs = [-1 1];
%     for i = 1:2
%         for j= 1:2
%             for k = 1:2
%                 u =  2^(-3/2) * exp(-1i * (pi / 4) * (signs(i) + signs(j) + signs(k))) ;
%                 Tr = Tr + u * det((I - Gs{i}) * (I - Gs{j}) * (I - Gs{k}) / 2^3 + ...
%                     expm(1i .* alpha .* N) * (I + Gs{i}) * (I + Gs{j}) * (I + Gs{k}) / 2^3);
%             end
%         end
%     end
    
    L2 = v2 - u2 + 1;
    r = Tr * exp(-1i * alpha * L2);
end