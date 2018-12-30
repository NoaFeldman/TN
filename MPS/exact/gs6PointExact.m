function gs6PointExact(intervals, gap, filename)
    cicjOrig = cicjForInfiniteEnv(gap * 2 + intervals(length(intervals))*3);
    alphas = -3.14:0.01:3.14;
    sFull = zeros(length(intervals), 1);
    s1Alpha = zeros(length(intervals), length(alphas));
    s2Alpha = zeros(length(intervals), length(alphas));
    s3Alpha = zeros(length(intervals), length(alphas));
    s4Alpha = zeros(length(intervals), length(alphas));
    s5Alpha = zeros(length(intervals), length(alphas));
    for i = 1:length(intervals)
        u1 = 1;
        v1 = u1 + intervals(i) - 1;
        u2 = v1 + gap + 1;
        v2 = u2 + intervals(i) - 1;
        u3 = v2 + gap +1;
        v3 = u3 + intervals(i) - 1;
        cicj = [cicjOrig(u1:v1, u1:v1) cicjOrig(u1:v1, u2:v2) cicjOrig(u1:v1, u3:v3); ...
                cicjOrig(u2:v2, u1:v1) cicjOrig(u2:v2, u2:v2) cicjOrig(u2:v2, u3:v3); ...
                cicjOrig(u3:v3, u1:v1) cicjOrig(u3:v3, u2:v2) cicjOrig(u3:v3, u3:v3)];
      
        [~, V] = eig(cicj);
        f = diag(V);
        s1Alpha(i, :) = getSAlpha(1, alphas, f, 1);
        s2Alpha(i, :) = getSAlpha(2, alphas, f, 1);
        s3Alpha(i, :) = getSAlpha(3, alphas, f, 1);
        s4Alpha(i, :) = getSAlpha(4, alphas, f, 1);
        s5Alpha(i, :) = getSAlpha(5, alphas, f, 1);
        sFull(i) = getExactEEForFlux(f, 0);
    end
    clear u1 v1 u2 v2 u3 v3 f V cicjOrig
    save(filename);
end