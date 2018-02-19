function A=skipzeros(A,varargin)
% function A=skipzeros(A [,eps])
% Wb,Nov10,07

  if nargin>1
     if nargin>2
        eval(['help ' mfilename]);
        if nargin || nargout, error('Wb:ERR','invalid usage'), end, return
     end
     eps=varargin{1};
     if ~isscalar(eps) || eps<0, error('Wb:ERR','invalid usage'); end
     if nargin>2
        eval(['help ' mfilename]);
        if nargin || nargout, error('Wb:ERR','invalid usage'), end, return
     end
     o={eps};
  else o={};
  end

% include possible set CG coefficients -- Wb,May14,10
  for k=1:numel(A)
   % QSpace() required by matlab/2016a // Wb,Aug03,16
     A(k)=QSpace(skipZerosQS(A(k),o{:}));
  end

end

