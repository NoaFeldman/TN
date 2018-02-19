function mvlabel(varargin)
% mvlabel - move label
%
% USAGE:
%
%    mvlabel ([ah,] lb, eps)
%
% where
%
%    ah  = axis handle
%    lb  = 'X|Y|Z|T'   - x label, y label, z label, title respectively
%    eps = shift in normalized units
%
% Wb,May22,03

  getopt('INIT',varargin);
     pos =getopt('pos',[]);
     POS =getopt('POS',[]);
  varargin=getopt('get_remaining');

  if ~isempty(POS)
   % set absolute position for all axes labels of current figure
     if ~ischar(varargin{1}), error('Wb:ERR',...
     'invalid usage (no axis handle with POS)'); end
     ah=findall(gcf,'type','axes'); n=numel(ah);
     for i=1:n
        if ~isempty(regexp(get(ah(i),'tag'),'legend|scribe')), continue; end
        mvlabel(ah(i),'pos',POS,varargin{:});
     end
   % NB! default vertical-alignment of labels: bottom
   % eg. for y-label, may include "_{ }" for proper alignment
   % with other axes sets
     return
  end
  
  narg=length(varargin);

  if isempty(pos)
     if narg==2, [ah,lb,eps]=deal(gca,varargin{1:2});
     else
        if narg~=3, eval(['help ' mfilename]); return; end
        [ah,lb,eps]=deal(varargin{1:3});
     end
  else
     if narg==1, ah=gca; lb=varargin{1};
     elseif narg==2, [ah,lb]=deal(varargin{1:2});
     else eval(['help ' mfilename]); return; end
     eps=pos;
  end

  ne=length(eps); unorm=1;
  if ne<1 || ne>3, eval(['help ' mfilename]); return; end

  switch lower(lb)
    case 'x',  h=get(ah,'xlabel'); if ne==1, eps=[0 eps 0]; end
    case 'y',  h=get(ah,'ylabel'); if ne==1, eps=[eps 0 0]; end
    case 'z',  h=get(ah,'zlabel'); if ne==1, eps=[eps 0 0]; end
    case 't',  h=get(ah,'title' ); if ne==1, eps=[0 eps 0]; end
  otherwise
    unorm=0;
    if isnumeric(lb) && isscalar(lb) && ...
       ishandle(lb) && isequal(get(lb,'type'),'text')

       h=lb; eps=[eps,0]; % expecting (dx,dy) data for eps
    else
       eval(['help ' mfilename]); error('Wb:ERR');
       return
    end
  end

  if iscell(h) % from multiple axis sets
  h=cat(1,h{:}); end

  if length(eps)==2, eps=[eps 0]; end

  if ~isempty(pos)
      set(h,'Units','norm','Pos',eps);
  else
      for i=1:numel(h)
         u=get(h(i),'Units'); % remember to reset below
         set(h(i),'Units','norm');
         set(h(i),'Pos',get(h(i),'Position')+eps);
         if ~unorm, set(h(i),'Units',u); end
      end
  end

end

