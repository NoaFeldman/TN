function h=poslabel(varargin)
% Function h=poslabel([ah,] n [,pos,opts])
%
%    position figure / panel label (a), ...
%
% Options
%
%    'dx',..  (in normalized units in current axis set)
%     remaining options are applied to text object.
%
% Wb,Apr09,11

  ah0=gca; 
  if nargin && isnumeric(varargin{1}) && isaxis(varargin{1})
       ah=reshape(varargin{1},1,[]); setax(ah(1)); varargin=varargin(2:end);
  else ah=ah0; end

  if numel(varargin)<1
     eval(['help ' mfilename]);
     if nargin || nargout, error('Wb:ERR','invalid usage'), end, return
  end

  getopt('init',varargin);
     rmflag=getopt('-rm');
     fix   =getopt('-fix');
     dx    =getopt('dx',[]);
  varargin=getopt('get_remaining'); narg=length(varargin);

  tag=mfilename; % 'poslabel'

  h=findall(ah,'type','text','tag',tag);
  if ~fix, delete(h); end

% always remove existing ones
  if ~narg
     if fix
        if numel(h)
           if numel(h)>1, error('Wb:ERR',...
            '\n   ERR got multiple labels !??'); end
           n=get(h,'string');
        else
           wblog('WRN','%s: no label to fix',mfilename);
           clear h, return
        end
     else
        if rmflag && ~nargout, clear h; return, end
        error('Wb:ERR','\n   ERR invalid usage');
     end
  else
     n=varargin{1};
     varargin=varargin(2:end); narg=numel(varargin);
  end

  p=[];
  if narg, q=varargin{1};
     if ischar(q) && numel(q) && numel(q)<=2 && isequal(q,upper(q))
     p=q; varargin=varargin(2:end); end
  end
  if isempty(p), p='NW'; end

  pos=postrans(p);
  if ~isempty(dx)
     if numel(dx)~=2, error('Wb:ERR','invalid dx'); end
     pos=pos+reshape(dx,1,[]);
  end

  if ischar(n)
     if ~isempty(n) && (n(1)=='(' || n(1)=='\' || n(1)==' ')
         s=n; else s=sprintf('(%s)',n); end
     h=postext(pos,s);
  else
     if ~isnumber(n), n, error('Wb:ERR','\n   ERR invalid usage'); end
     h=postext(pos,sprintf('(%s)',char('a'+n-1)));
  end

  set(h,'tag',tag,'BackgroundC','w','Margin',0.1,'FontSize',16,varargin{:});
  if ~nargout, clear h; end

  if ~isequal(ah,ah0), setax(ah0); end

end

