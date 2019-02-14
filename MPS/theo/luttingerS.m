function S = luttingerS(N, cici, cicicjcj)
    S = 0;
    q = 2*pi./N;
    for i = 1:N
        for j = 1:N
            S = S + 1./N .* exp(1i .* q * (i - j)) * ...
                (cicicjcj(i, j)  - cici(i)*cici(j));
        end
    end
end