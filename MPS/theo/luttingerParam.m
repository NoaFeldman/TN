function K = luttingerParam(N, cici, cicicjcj)
    S = luttingerS(N, cici, cicicjcj);
    q = 2*pi./N;
    K = pi * S / q;
end