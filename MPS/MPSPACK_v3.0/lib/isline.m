function r=isline(hh)
% Functions: i=isline(h)
% 
%    Check whether given input is valid line handles
% 
% Wb,Nov12,06

  if nargin~=1, eval(['help ' mfilename]); return; end
  if ~ishandle(hh), r=0; return; end

  m=numel(hh); r=1;

  for i=1:m
     if ~ishandle(hh(i)) || ~isequal(get(hh(i),'Type'),'line') && ...
        ~isequal(get(hh(i),'Type'),'hggroup') % errorbar
     r=0; return; end
  end

end

