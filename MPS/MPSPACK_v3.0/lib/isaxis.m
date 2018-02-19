function ia=isaxis(hh)
% Functions: isaxis(h)
% 
%    Check whether given input represents valid axis handles
% 
% Wb,Nov12,06

  if nargin~=1, eval(['help ' mfilename]); return; end
  if any(ishandle(hh)) && ~isempty(hh), ia=1; else ia=0; return; end

  s=size(hh); m=prod(s);
  for i=1:m
     if ~ishandle(hh(i)) || ~isequal(get(hh(i),'Type'),'axes')
     ia=0; return; end
  end

end

