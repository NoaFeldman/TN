function saveR3Negativity(qs, L, l1Overl2, filename)
    l2 = L / (1 + l1Overl2);
    l1 = L - l2;
    r3 = zeros(1, length(qs));
    cicj = getCiCj0Matrix(1000);
    for i = 1:L
        for j = 1:L
            if (i == j)
                cicj(i, j) = 1/2;
            else
                cicj(i, j) = sin(pi * (i-j) / 2)/(pi * (i-j));
            end
        end
    end
    for i = 1:length(qs)
        r3(i) = chargeRenyi3(cicj, qs(i), 1, l1, L);
    end
    save(strcat(filename, '_', int2str(l1Overl2)), 'r3', 'qs');
end
        