function v=getuser(varargin)
% Function u=getuser([hh,] tag)
%
%    Get field tag of userdata of current axes (or hh).
%
% Wb,Sep02,08

  if nargin && ishandle(varargin{1})
     ah=varargin{1}; varargin=varargin(2:end);
     narg=length(varargin);
  else
     ah=gca; narg=nargin;
  end

  getopt('init',varargin);
     rflag=getopt('-rm');
  varargin=getopt('get_remaining'); narg=length(varargin);

  if narg==0 && numel(ah)==1
     v=get(ah,'UserD');
     return
  end

  if narg~=1 || ~ischar(varargin{1})
     eval(['help ' mfilename]);
     if nargin || nargout, error('Wb:ERR','invalid usage'), end, return
  end

  f=varargin{1};
  ah=reshape(ah,[],1); n=length(ah);

  for i=1:n, h=ah(i);
     u=get(h,'UserData');
     if ~isempty(h) && isfield(u,f)
        v=getfield(u,f);
        if rflag, set(h,'UserD',rmfield(u,f)); end
     else v=[]; end
  end

end

