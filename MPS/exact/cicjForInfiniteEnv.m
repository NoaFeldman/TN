function cicj = cicjForInfiniteEnv(L)
    cicj = zeros(L);
    for i = 1:L
        cicj(i, i) = 1/2;
        for j = i + 1:L
            cicj(i, j) = sin(pi * (i-j) / 2)/(pi * (i-j));
            cicj(j, i) = cicj(i, j);
        end
    end
end
    