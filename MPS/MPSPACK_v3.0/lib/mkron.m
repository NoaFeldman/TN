function XX = mkron(varargin)
% Function: XX = mkron(varargin)
%
%    Extends matlabs kron() to more than one two arguments
%    NB! first index is assumed fastest.
%
% Options (last argument)
%
%   'rowmajor'  take last index fastest
%
% Wb,Jun16,07
   
  if nargin<2
     if nargin && isnumeric(varargin{1}), XX=varargin{1}; return; end
     eval(['help ' mfilename]);
     if nargin || nargout, error('Wb:ERR','invalid usage'), end, return
  end

  rm=0;
  if isequal(varargin{end},'rowmajor')
     rm=1; varargin=varargin(1:end-1);
  end

  n=length(varargin);
  for i=1:n
     if ~isnumeric(varargin{i})
     error('Wb:ERR','input must be numeric'), end
  end

  if rm
     XX=varargin{1};
     for i=2:n, XX=kron(XX,varargin{i}); end
  else
     XX=varargin{n};
     for i=n-1:-1:1, XX=kron(XX,varargin{i}); end
  end

end

