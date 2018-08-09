function cicj = getCiCj0Matrix(L)
    ckcq = zeros(L, L);
    for i = 1 : L/2
        ckcq(i, i) = 1;
    end
    S = realSpaceToDualSpace(L); 
    cicj = S * ckcq * S';
end