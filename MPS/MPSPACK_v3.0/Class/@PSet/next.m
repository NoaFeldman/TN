function i=next(P)
% Function next(P,i)
% Wb,Aug04,08

  if nargin>1
     eval(['help ' mfilename]);
     if nargin || nargout, error('Wb:ERR','invalid usage'), end, return
  end

  if P.i==P.n, i=0; else P.i=P.i+1; i=P.i; end
   
  assignin('caller',inputname(1),P);

end

