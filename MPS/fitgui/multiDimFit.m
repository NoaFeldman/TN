function [a,aerr,cov,chisq,yfit] = multiDimFit(funcName, x, y, aGlobal, aVariant, fixedGlobal, fixedVariant)
    % Fit func to x with multiple y sets, aGlobal are parameters global to
    % all y sets, and aVariant varying by y set.
    % Assuming y is a nXd matrix, with n = x size and d = dimension
    % num.
    % aGlobal in func must be prior to aVariant in a.
    % For now, fixedVariant only handles 1 variable!!!!!!!!!!!!!!!!!
    
    % Assuming x is a row vector
    n = length(x);
    d = length(y(1, :));
    multiX = repmat(x, [1, d]);
    multiY = reshape(y, [1, n*d]);
    a = [aGlobal, repmat(aVariant, [1, d])];
    
    fid = fopen('multiFunc.m', 'wt');
    fprintf(fid, 'function res = multiFunc(x, a, fixed, fixedVariant)\n');
    fprintf(fid, 'res = zeros(length(x), 1);\n');
    fprintf(fid, 'for i = 0:%d\n', d-1);
    fprintf(fid, 'res(i * %d + 1:(i+1) * %d) = %s(x(i * %d + 1:(i+1) * %d), [a(1:%d)         a(%d   + 1 + i * %d              :%d   + (i+1) *   %d)], [fixed fixedVariant(i+1)]);\n', ...
                           n,            n,funcName,     n,             n,  length(aGlobal), length(aGlobal), length(aVariant),length(aGlobal), length(aVariant));
    fprintf(fid, 'end');
    
    [a,aerr,cov,chisq,yfit] = fitnonlin(multiX, multiX, multiY, multiX*0.01, multiY*0.01, 'multiFunc', a, fixedGlobal, fixedVariant);
    delete multiFunc.m
end
          