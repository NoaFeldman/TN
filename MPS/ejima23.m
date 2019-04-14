Ns = 32:32:1024;
for i = 1:length(Ns)
    N = Ns(i);
    H = zeros(N);
    for j = 1:N-1
        H(j, j+1) = 1;
        H(j+1, j) = 1;
    end
    [vecs, vals] = eig(H);
    vals = diag(vals);
    E0O(i) = sum(vals(1:N/2));
    E1O(i) = sum(vals(1:N/2 + 1));
    
    H(1, N) = 1;
    H(N, 1) = 1;
    [vecs, vals] = eig(H);
    
    vals = diag(vals);
    E0(i) = sum(vals(1:N/2));
    E2(i) = sum(vals(1:N/2 + 2));
end
    