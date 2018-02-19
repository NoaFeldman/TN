function bt(n)
% Function bt() - backtrace (when debugging)
% Wb,Aug27,08
   
  if ~nargin, n=0; end
% [S,i]=dbstack; will always take *this (=bt.m) as current stack line => i=1
  [S]=dbstack('-completenames'); S=S(max(1,2+n):end); if isempty(S), return; end

% -------------------------------------------------------------------- %
% get current frame
% => useless, since for the call to dbstack here will always
%    return *this as current frame!
% -------------------------------------------------------------------- %
%    [q]=evalc('dbstack'); q=[10 q(1:end-1)];
%     i=find(q==10); i=i(max(1,2+n):end);
%   % [q(i+1); q(i+2); q(i+3);; q(i+4); q(i+5)]'
%     j=find(q(i+1)=='>');
%     if ~isempty(j), S(j(1)).curr=1; end
% -------------------------------------------------------------------- %

  dispstack(S)

end

