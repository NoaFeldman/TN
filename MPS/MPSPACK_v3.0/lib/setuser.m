function u=setuser(varargin)
% Function setuser([hh,] tag, value)
%
%    set field <tag> of userdata of given handles (default: current axis)
%    to value specified.
%
% Wb,Sep02,08

  if nargin && ishandle(varargin{1})
     ah=varargin{1}; varargin=varargin(2:end);
     narg=length(varargin);
  else
     ah=gca; narg=nargin;
  end

  if narg~=2
     eval(['help ' mfilename]);
     if nargin || nargout, error('Wb:ERR','invalid usage'), end, return
  end

  f=varargin{1}; v=varargin{2};
  ah=reshape(ah,[],1); n=length(ah);

  for i=1:n, h=ah(i);
     u=get(h,'UserData');
     if isstruct(u), u=setfield(u,f,v); else u=struct(f,{v}); end
     set(h,'UserData',u);
  end

  if ~nargout, clear u; end

end

