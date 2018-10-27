function [a,aerr,cov,chisq,yfit] = multiDimFit(funcName, x, y, aGlobal, aVariant, fixed)
    % Fit func to x with multiple y sets, aGlobal are parameters global to
    % all y sets, and aVariant varying by y set.
    % Assuming y is a nXd matrix, with n = x size and d = dimension
    % num.
    % aGlobal in func must be prior to aVariant in a.
    
    % Assuming x is a row vector
    n = length(x);
    d = length(y(1, :));
    multiX = repmat(x, [1, d]);
    multiY = reshape(y, [1, n*d]);
    a = [aGlobal, repmat(aVariant, [1, d])];
    
    fid = fopen('multiFunc.m', 'wt');
    fprintf(fid, 'function res = multiFunc(x, a, fixed)\n', 'tst');
    fprintf(fid, 'res = zeros(length(x), 1);');
    fprintf(fid, 'for i = 0:%d\n', d-1);
    fprintf(fid, 'res(i * %d + 1:(i+1) * %d) = %s(x(i * %d + 1:(i+1) * %d), [a(1:%d)         a(%d   + 1 + i * %d              :%d   + (i+1) *   %d)], fixed);\n', ...
                           n,            n,funcName,     n,             n,  length(aGlobal), length(aGlobal), length(aVariant),length(aGlobal), length(aVariant));
    fprintf(fid, 'end');
    
    [a,aerr,cov,chisq,yfit] = fitnonlin(multiX, multiX, multiY, multiX*0.01, multiY*0.01, 'multiFunc', a, fixed);
    delete multiFunc.m
end
          