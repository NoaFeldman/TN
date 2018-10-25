clear r3;
L = 10000;
ls = 500:100:1000;
for i = 1:length(ls)
    l = ls(i);
    r3(i, :) = exactNegAfterQuenchFlux(L, L/2 - l + 1, L/2, L/2 + l, 0, l, 1, pwd, 'tst', 4).';
end
save('negsT0TLong', 'ls', 'r3', 'L');
