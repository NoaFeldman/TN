function printf(varargin)
% printf - emulating regular printf
%
%    Usage : printf(fmt, ...)
%    Equivalent to disp(sprintf(...))
%
% Wb,Aug13,05

% disp(sprintf(varargin{:}));
  fprintf(1,varargin{:});
  fprintf(1,'\n');

end

