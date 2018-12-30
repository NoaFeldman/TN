function sn = gs6PointTheo(alphas, a, fixed)
    s = size(alphas);
    if s(1) == 1
        alphas = alphas.';
    end
    n = fixed(1);
    u1 = 1;
    interval = fixed(2);
    gap = fixed(3);
    v1 = u1 + interval - 1;
    u2 = v1 + gap + 1;
    v2 = u2 + interval - 1;
    u3 = v2 + gap +1;
    v3 = u3 + interval - 1;
    usAndVs = [u1 v1 u2 v2 u3 v3];
    nuc = abs(a(1));
    K = 1;
    order1 = ones(length(alphas), 1);
    for i = 1:length(usAndVs)
        signI = (mod(i, 2) - 1/2)*2;
        for j = i+1:length(usAndVs)
            signJ = (mod(j, 2) - 1/2)*2;
            order1 = order1 .* (nuc*(usAndVs(j) - usAndVs(i))).^(K/n .* signI .* signJ .* (alphas./(2*pi)).^2);
        end
    end
    sn = order1;
end