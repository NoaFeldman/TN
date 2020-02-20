Ns = 256:256:4096;
for i = 1:length(Ns)
    N = Ns(i);
    H = zeros(N);
    for j = 1:N-1
        H(j, j+1) = 1;
        H(j+1, j) = 1;
    end
    [vecs, vals] = eig(H);
    vals = diag(vals);
    E0Open(i) = sum(vals(find(vals < 0)));
end

for i = 1:length(Ns)
    N = Ns(i);
    H = zeros(N);
    for j = 1:N-1
        H(j, j+1) = 1;
        H(j+1, j) = 1;
    end
    H(1, N) = 1;
    H(N, 1) = 1;
    [vecs, vals] = eig(H);
    vals = diag(vals);
    E0Periodic(i) = sum(vals(find(vals < 0)));
end