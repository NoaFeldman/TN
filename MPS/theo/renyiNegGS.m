function res = renyiNegGS(alpha, a, fixed)
    nuc = a(1);
    w = a(2);
    l1 = fixed(1);
    l2 = fixed(2);
    n = fixed(3);
    K = fixed(4);
    
    l1 = nuc*l1;
    l2 = nuc*l2;
    alpha1 = alpha / (2*pi);
    alpha2 = -2 .* alpha / (2*pi);
    alpha3 = alpha / (2*pi);
    
    order1 = l1.^(2 * K / n .* alpha1 .* alpha2) .* l2.^(2 * K / n .* alpha2 .* alpha3) .* (l1 + l2).^(2 * K / n .* alpha1 .* alpha3);
    order2 = 0;
    order2Shift = [0, 1, -1];
    for i = 1:3
        for j = 1:3
            for k = 1:3
                if i ~= j && i ~=k && k ~= j
                    order2 = order2 + ...
                        w .* l1.^(2 * K / n .* (alpha1 + order2Shift(i)) .* (alpha2 + order2Shift(j))) .* ...
                        l2.^(2 * K / n .* (alpha2 + order2Shift(j)) .* (alpha3 + order2Shift(k))) .* ...
                        (l1 + l2).^(2 * K / n .* (alpha1 + order2Shift(i)) .* (alpha3 + order2Shift(k)));
                end
            end
        end
    end
    res = order1 + order2;
end