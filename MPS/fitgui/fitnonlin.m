function [a,aerr,cov,chisq,yfit] = fitnonlin(x, x_res, y, sigx, sigy, fitfun, a0, fixed, fixedVariant)

% Version 2 (creatd by Adiel Meyer 15.10.2012)
%
% FITNONLIN Fit a nonlinear function to data.
%    [a,aerr,chisq,yfit] = fitnonlin(x,y,sig,fitfun,a0) 
%
%    Inputs:  x -- the x data to fit
%             xres -- a higher resolution version of x for plotting, where
%             yfit is plotted.
%             y -- the y data to fit
%             sigx -- the uncertainties on the data points
%             sigy -- the uncertainties on the data points
%             fitfun -- the name of the function to fit to
%             a0 -- the initial guess at the parameters 
%
%    Outputs: a -- the best fit parameters
%             aerr -- the errors on these parameters
%             chisq -- the final value of chi-squared
%             yfit -- the value of the fitted function
%                     at the points in x_res
%  
%    Note: "fitfun" should be in a .m file similar to the
%    following example.
%
%          The following lines are saved in a file called
%          "sinfit.m", and the routine is invoked with
%          the fitfun parameter equal to 'sinfit' (including
%          the quotes)
%
%          function f = sinfit(x,a)
%          f = a(1)*sin(a(2)*x+a(3));
%

% first set up the parameters needed by the algorithm
x = x(:);
y = y(:);
sigx = sigx(:);
sigy = sigy(:);


stepdown = 0.1;
stepsize = abs(a0)*0.01+eps ; % the amount each parameter will be varied by in each iteration
chicut = 0.00001;  % maximum differential allowed between successive chi sqr values

% These parameters can be varied if you have reason to believe your fit is
% converging to quickly or that you are in a local minima of the chi
% square.

a = a0;
fprintf(1, 'Iteration No.       ');
iter=0;
%chi2 = calcchi2(x,y,sig,fitfun,a);
chi2 = calcchi2(x,y,sigx,sigy,fitfun,a,fixed, fixedVariant);
chi1 = chi2*2;

% keep looking while the value of chi^2 is changing

while (abs(chi2-chi1))>chicut

  [anew,stepsum,stopflag,iter] = gradstep(x,y,sigx,sigy,fitfun,a,stepsize,stepdown,iter,fixed, fixedVariant);
  a = anew;
  stepdown = stepsum;
  chi1 = chi2;
  chi2 = calcchi2(x,y,sigx,sigy,fitfun,a,fixed, fixedVariant);
  if stopflag==1
	fprintf(2,'\n can''t minimize chi^2 \n try different initial parameters\n')
	break, end
end

% calculate the returned values
[aerr,cov] = sigparab(x,y,sigx,sigy,fitfun,a,stepsize,fixed, fixedVariant);
chisq = calcchi2(x,y,sigx,sigy,fitfun,a,fixed, fixedVariant);
if (isempty(fixedVariant))
    yfit = feval(fitfun,x_res,a,fixed);
else
    yfit = feval(fitfun,x_res,a,fixed, fixedVariant);
end

%----------------------------------------------------------------------- 
% the following function calculates the (negative) chi^2 gradient at
% the current point in parameter space, and moves in that direction
% until a minimum is found
% returns the new value of the parameters and the total length travelled

function [anew,stepsum,stopflag,iter] = gradstep(x,y,sigx,sigy,fitfun,a,stepsize, stepdown,iter,fixed, fixedVariant)
stopflag=0;

chi2 = calcchi2(x,y,sigx,sigy,fitfun,a,fixed, fixedVariant);

grad = calcgrad(x,y,sigx,sigy,fitfun,a,stepsize,fixed, fixedVariant);

chi3 = chi2*1.1;

chi1 = chi3;

% cut down the step size until a single step yields a decrease 
% in chi^2
stepdown = stepdown*2;
j = 0;
maxiter=100000;
if (chi3 < chi2)
    k = 1;
end
while chi3>chi2
  stepdown = stepdown/2;
  anew = a+stepdown*grad;
  chi3 = calcchi2(x,y,sigx,sigy,fitfun,anew,fixed, fixedVariant);
  j=j+1;
  iter=iter+1;
  fprintf('\b\b\b\b\b\b');
  fprintf(1, '%6d', iter);
  
  if (iter > maxiter) 
      stopflag=1;
      break;
  end;
end

stepsum = 0;

% keep going until a minimum is passed

while chi3<chi2
%   fprintf(['\n ',num2str(anew),' \n'])
  stepsum = stepsum+stepdown;
  chi1 = chi2;
  chi2 = chi3;
  anew = anew+stepdown*grad;
  chi3 = calcchi2(x,y,sigx,sigy,fitfun,anew,fixed, fixedVariant);
  iter=iter+1;
  fprintf('\b\b\b\b\b\b');
  fprintf(1, '%6d', iter);
  
  if (iter > maxiter) 
      stopflag=1;
      break; 
  end;
  %fprintf(2, 'iteration 2 No. %d\n', t)
  %fprintf(2,'can not find fit parameters, try differrent intial parameters\n')
end

% approximate the minimum as a parabola
step1 = stepdown*((chi3-chi2)/(chi1-2*chi2+chi3)+.5);
anew = anew - step1*grad;


%--------------------------------------------------------------
% this function calculates the (negative) gradient at a point in 
% parameter space

function grad = calcgrad(x,y,sigx,sigy,fitfun,a, stepsize, fixed, fixedVariant)

f = 0.01;
[dum, nparm] = size(a);

grad = a;
chisq2 = calcchi2(x,y,sigx,sigy,fitfun,a, fixed, fixedVariant);  
for i=1:nparm

  a2 = a;
  da = f*stepsize(i);
  a2(i) = a2(i)+da;
  chisq1 = calcchi2(x,y,sigx,sigy,fitfun,a2, fixed, fixedVariant);
  grad(i) = chisq2-chisq1;

end

t = sum(grad.^2);
grad = stepsize.*grad/sqrt(t);

%------------------------------------------------------------
% this function calculates the errors on the final fitted 
% parameters by approximating the minimum as parabolic
% in each parameter.
function [err1,cov]=sigparab(x,y,sigx,sigy,fitfun,a,stepsize,fixed, fixedVariant)

[dum, nparm] = size(a);

for j=1:nparm

    da(j) = stepsize(j);

    a1=a;
    a1(j)=a1(j)+da(j);
   
     a2= a;
     a3= a1;
     for k=1:nparm
          da(k) = stepsize(k);
         a2(k)=a2(k)+da(k);
         a3(k)=a3(k)+da(k); 
          dCHI2da(j,k)=0.5*(calcchi2(x,y,sigx,sigy,fitfun,a,fixed, fixedVariant)-calcchi2(x,y,sigx,sigy,fitfun,a1,fixed, fixedVariant)-calcchi2(x,y,sigx,sigy,fitfun,a2,fixed, fixedVariant)+calcchi2(x,y,sigx,sigy,fitfun,a3,fixed, fixedVariant))/da(j)/da(k);
     end
 end
errc=inv(dCHI2da);
cov=errc;
for j=1:nparm
     err1(j) = sqrt(abs(errc(j,j)));
end

  