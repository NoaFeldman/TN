%------------------------------------------------------------
% this function just calculates the value of chi^2
function chi2 = calcchi2(x,y,sigx,sigy,fitfun,a, fixed, fixedVariant)
%chi2 = sum( ((y-feval(fitfun,x,a)) ./sig).^2);
xup = x + sigx;
xdwn = x - sigx;
chi2 = sum( ((y-feval(fitfun,x,a, fixed, fixedVariant)).^2)./ ...
    (sigy.^2 + ((feval(fitfun,xup,a, fixed, fixedVariant) - feval(fitfun,xdwn,a, fixed, fixedVariant))./2).^2) );
if chi2 < 0
    k = 1;
end