clear rQ0 rQ1 rQ2 rQ3 rQ4 rQ5 rQ6 rQ7 rQ8 rQ9 rAlpha
% dirName = 'negsL256Delta-05';
% tSteps = 0:2:44;
dirName = 'quenchNegsl64Delta-03';
tSteps = 0:2:128;
t = tSteps*0.25;
alphas = -3.14:0.01:3.14;
zeroIndex = find(alphas == 0);
rAlpha = zeros(length(t), length(alphas));
traceNorm = zeros(length(t), 10);
% n = 1;
for i = 1:length(tSteps)
    tStep = tSteps(i);
    load(strcat(dirName, '/step', int2str(tStep)));
    temp = arrangeNegativityDMRGRenyi(res.l64, n);
    rQ0(i) = temp.rnq(1); % ./ temp.rnalpha(zeroIndex);    
    rQ1(i) = temp.rnq(2); % ./ temp.rnalpha(zeroIndex);    
    rQ2(i) = temp.rnq(3); % ./ temp.rnalpha(zeroIndex);    
    rQ3(i) = temp.rnq(4); % ./ temp.rnalpha(zeroIndex);
    rQ4(i) = temp.rnq(5); % ./ temp.rnalpha(zeroIndex);
    rQ5(i) = temp.rnq(6); % ./ temp.rnalpha(zeroIndex);    
    rQ6(i) = temp.rnq(7); % ./ temp.rnalpha(zeroIndex);
    rQ7(i) = temp.rnq(8); % ./ temp.rnalpha(zeroIndex);    
    rQ8(i) = temp.rnq(9); % ./ temp.rnalpha(zeroIndex);
    rQ9(i) = temp.rnq(10);% ./ temp.rnalpha(zeroIndex);
    rAlpha(i, :) = temp.rnalpha;
    traceNorm(i, :) = temp.traceNorm(1:10);
    truncErrs(i) = res.truncErrs;
end