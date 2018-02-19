function A=diag(A,varargin)
% Function: A=diag(A [,OPTS])
%
%    extract diagonal of operator in QSpace format
%    alternatively, if only diagonal is stored,
%    diagonals become full again.
%
% Options
%
% Wb,Sep08,06

  getopt('init',varargin);
     dflag=getopt('-d'); % full data vector only
     trans=getopt('-t'); % transpose flag
  getopt('check_error');

  if nargin<1 || nargout>1
     eval(['help ' mfilename]);
     if nargin || nargout, error('Wb:ERR','invalid usage'), end, return
  end

  if ~mpsIsQSpace(A) || prod(size(A))~=1
  error('Wb:ERR','Invalid QSpace (need single QSpace object)'); end

  if isempty(A.Q)
     if isempty(A.data)
        if dflag, A=[]; end
        % else just return empty QSpace - ok.
     else
        if numel(A.data)~=1
           error('Wb:ERR','\n   ERR unexpected QSpace structure'); end
        q=diag(A.data{1}); if trans, q=q.'; end
        if dflag, A=q; else A.data{1}=q; end
     end
     return
  end

  if length(A.Q)~=2, wblog('ERR',...
     '%s requires rank-2 object (%d).',mfilename,length(A.Q)); end
  if ~isequal(A.Q{:}) && ~isequal(A.Q{1},fliplr(A.Q{2})), wblog('ERR',...
     '%s requires block-diagonal operator.',mfilename); end

  n=length(A.data);
  if trans
       for i=1:n, A.data{i}=diag(A.data{i}).'; end
  else for i=1:n, A.data{i}=diag(A.data{i})  ; end
  end

  if dflag
     if trans, i=2; else i=1; end
     A=cat(i,A.data{:});
  end

end

