function P=reset(P,i)
% Function reset(P,i)
% Wb,Aug04,08

  if nargin>2 || nargin==2 && ~isscalar(i)
     eval(['help ' mfilename]);
     if nargin || nargout, error('Wb:ERR','invalid usage'), end, return
  end

  if nargin<2, i=1; end
  P.i=i;
   
  if nargout==0, assignin('caller',inputname(1),P); clear P; end

end

