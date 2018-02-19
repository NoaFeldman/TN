function [i,bfac,e]=sameup2fac(A,B,varargin)
% function [i,bfac,e]=sameup2fac(A,B [,eps])
%
%   checks whether A = bfac * B within numerical noise eps (1E-12).
%
% Options:
%
%   '-ev'  check whether  A = B * Lambda, i.e. whether A and B
%          are related like eigenvectors.
% 
% Wb,Jun29,10

  if nargin<2
     eval(['help ' mfilename]);
     if nargin || nargout, error('Wb:ERR','invalid usage'), end, return
  end

  getopt('init',varargin);
     evflag=getopt('-ev');
  eps=getopt('get_last',1E-10);

  if ~isnumeric(A) || ~isnumeric(B)
  error('Wb:ERR','invalid input (numeric data expected)'); end

  i=0; e=nan;
  if ~isequal(size(A),size(B)), return; end

  if evflag
   % A = B Lambda, while assuming that columns in B are orthogonal
     B2=B'*B; Lambda=inv(B2)*(B'*A);
     e=norm(Lambda-diag(diag(Lambda)))/norm(Lambda);
     if e>eps, return; end; e=norm(A-B*Lambda)/norm(B2);
     if e>eps, return; end
  else
     a=A(:); b=B(:); % |a> = |b><b|a> / <b|b>  == bfac |b>

     b2=b'*b; bfac=full(real(b'*a)/b2);

     e=norm(a-bfac*b)/sqrt(b2);
     if e<1E-8 && e>1E-12, wblog('WRN','got small deviation (%.3g)',e); end
     if e>eps, return; end
  end

  i=1;

end

