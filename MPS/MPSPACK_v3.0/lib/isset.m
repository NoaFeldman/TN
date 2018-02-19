function i=isset(var,i0)
% function isset(var [,i0])
%
%    check whether variable 'var' is defined in current workspace.
%    if variable is not found, then use value of i0, instead (if specified).
%
% Wb,Feb03,10

% tags: gotq exists

  if nargin<1 || nargin>2 || ~ischar(var)
     helpthis, if nargin || nargout
     error('Wb:ERR','invalid usage'), end, return
  end

  tag='gotq__'; setuser(groot,tag,tag);

  evalin('caller',...
    ['if exist(''' var ''',''var''), setuser(groot,''' tag ''', ' var '); end']);
  q=getuser(groot,tag,'-rm');

  if ~isequal(q,tag) % variable exists
   % NB! isequal(48,'0') yields 1 (!) // Wb,Sep04,15
     if isempty(q ) || isequal(q, 0) || ischar(q) && isequal(q, '0'), i=0; else i=1; end
  elseif nargin>1 % variable does not exist: check default value
     if isempty(i0) || isequal(i0,0) || ischar(q) && isequal(i0,'0'), i=0; else i=1; end
  else i=0; end

end

