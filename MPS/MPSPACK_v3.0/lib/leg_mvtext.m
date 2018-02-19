function hh=leg_mvtext(l0,dy)
% Function: leg_mvtext(l0,dy)
%
%    with dy either a scalar (for all legend entries
%    or a vector for every legend entry
%
% Wb,May07,07

  if nargin~=2 || numel(l0)~=1 || ~ishandle(l0)
     helpthis, if nargin || nargout
     error('Wb:ERR','invalid usage'), end, return
  end

  t=get(l0,'type');
  if isequal(t,'axes')
     hh=findall(l0,'type','text','tag','')';
  elseif isequal(t,'legend')
     hh=get(l0,'ItemText'); % Wb,Aug03,16
  else l0, t, error('Wb:ERR','\n   ERR invalid handle'); 
  end

  m=numel(hh);

  if length(dy)==1, dy=dy(ones(m,1));
  elseif length(dy)>length(hh), error('Wb:ERR',...
  'invalid usage (%g/%g)',length(dy),length(hh)); end

  n=numel(dy);

  for i=1:n, h=hh(end-i+1);
     set(h,'Pos',get(h,'Pos')+[0 dy(i) 0]);
     setuser(h,'leg_dy',dy(i));
  end

  if ~nargout, clear hh; end

end

