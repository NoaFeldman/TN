function A=itagrep(A,varargin)
% function A=itagrep(A [,idx], <arguments to regepxrep>)
%
%    apply given regexprep() pattern to specified itags
%    if idx is not specified, then the itags for all dimensions
%    are considered.
%
% Wb,Jan14,15

  if nargin<3 || numel(A)~=1
     helpthis, if nargin || nargout
     error('Wb:ERR','invalid usage'), end, return
  end

  t=getitags(A);

  if isnumeric(varargin{1})
       it=varargin{1}; varargin=varargin(2:end);
  else it=1:numel(t); end

  t(it)=regexprep(t(it),varargin{:});
  A.info.itags=t; % clear t

  if ~nargout, vname=inputname(1);
     if ~isempty(vname)
        assignin('caller',vname,A);
        clear A
     end
  end

end

