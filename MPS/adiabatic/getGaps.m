hi=0;
JZi = 1;
JPMi = 0;
hf = 0;
JZf = 0;
JPMf = 1;
stepNum = 100;
N = 100;
minGap = 100;
for step = 0:100
    s = step/stepNum;
    [psi0, ~, ~, ~, E0] = getGroundState(N, hi * (1-s) + hf * s, + ...
        JPMi * (1-s) + JPMf * s, JZi * (1-s) + JZf * s, 0, 0, 0, 'open');
    % Start from 
    psiStart = getOrthogonalState(psi0);
    [psi1, ~, ~, ~, E1] = getGroundState(N, hi * (1-s) + hf * s, + ...
        JPMi * (1-s) + JPMf * s, JZi * (1-s) + JZf * s, 0, 0, 0, 'open', psiStart);
    gap = E1 - E0;
    if gap < minGap
        minGap = gap;
    end
end

