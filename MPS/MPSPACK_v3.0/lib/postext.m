function h=postext(varargin)
% Function h=postext(nx,ny, <fmt, vars, ... forwarded to printf> [,{topts}])
%
%    puts 'text' in the current axis at position (nx,ny)
%    in normalized units with (0,0) the lower left corner
%
%    Instead of nx,ny also a position similaro to the legend's
%    location flag can be specified: 'North','NorthEast',...
%
%    <fmt, vars, ...> is forwarded to printf
%
%    the last argument, if it is a cell, is interpreted as
%    text options applied to the text handle.
%
% Wb 2000, Wb,Mar31,07

  tag='postext';
  if nargin==1 && isequal(varargin{1},'-rm')
     delete(findall(gca,'tag',tag));
     return
  end

  if nargin<2, eval(['help ' mfilename]); return; end

  o={};

  if isnumeric(varargin{1})
     if isnumeric(varargin{2})
          x=varargin{1};    y=varargin{2};    varargin=varargin(3:end);
     else x=varargin{1}(1); y=varargin{1}(2); varargin=varargin(2:end);
     end
  else
     [pos,o]=postrans(varargin{1}); x=pos(1); y=pos(2);
     varargin=varargin(2:end);
  end

  getopt('init',varargin);
     rmflag=getopt('-rm');
     dx    =getopt('dx',0);
     totop =getopt('-top');
  varargin=getopt('get_remaining'); narg=length(varargin);

  if isempty(varargin), eval(['help ' mfilename]); return; end
  if ~isempty(dx)
     if length(dx)>2, error('Wb:ERR','invalid dx'); end
     x=x+dx(1); if length(dx)>1, y=y+dx(2);
  end

  if rmflag
     delete(findall(gca,'type','text','tag',tag));
     if ~narg, return; end
  end

  topts={};
  if iscell(varargin{end})
     topts=varargin{end}; varargin=varargin(1:end-1);
  end

  if length(varargin)>1
       fmt=regexprep(varargin{1},'\\+','\\\\');
       fmt=regexprep(fmt,'\\\\n','\\n'); % keep newlines!
       s=sprintf(fmt,varargin{2:end});
  elseif iscell(varargin{1})
     % concatenate cell of strings separated by newline
       s=strhcat(varargin{1},char(10));
  else s=varargin{1}; end
  s=regexprep(s,'\\\\+','\\');

  h=text(x,y,s,'Units','normalized','tag',tag,o{:},topts{:});

  if totop
     set(h,'Units','Data'); p=get(h,'Pos');
     p(3)=0.1; set(h,'Pos',p);
   % units=normalized, lets text vanish underneath patches for example
   % set(h,'Units','normalized');
  end

  if isequal(get(h,'HorizontalAl'),'right')
  set(h,'String',strvcat2(get(h,'String'))); end

  if ~nargout, clear h; end

end

