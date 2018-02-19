function dd=chopd(dd,n)
% function dd=chopd(dd [,n])
%
%    chop noise from double precision data
%    by converting to single and back to double.
%
%    if second argument is specified, round to n digits.
%
% See also chop.m
% Wb,Nov05,08

% NB! for complex numbers, rather use chop() than *this.

 % a=2*max(abs(dd(:))); % take twice to have real offset
 % add2 numerical noise!
   a=1; if ~isreal(dd), a=complex(a,a); end

   if nargin<2
      i=find(dd>0); if ~isempty(i), dd(i)=double(single(dd(i)+a))-a; end
      i=find(dd<0); if ~isempty(i), dd(i)=double(single(dd(i)-a))+a; end
   else
      fac=10^n;
      dd=round(dd*fac)/fac;
   end

end

