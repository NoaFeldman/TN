function lh=blurl(varargin)
% function blurl([lh, opts])
%
%    blur line by making it brighter with linewidth 2
%
% Options
% 
%    lh     set of line handles (if not provide, then all lines
%           within current axis set are taken)
% 
%   'cfac'  factor of how much to blur (0=no blur; 1=white) (0.8)
%   'cmin'  min of color to keep (relative to sqrt(white)=sqrt(3))
%   'linew' line width (2)
%   'undo'
%
% NB !any remaining options are directly applied to linehandles.
% 
% Wb,Oct31,06

  if ~nargin
     helpthis, if nargin || nargout
     error('Wb:ERR','invalid usage'), end, return
  end
  if isempty(varargin{1}), return; end

  bflag=0;

  if ishandle(varargin{1})
       lh=varargin{1}; varargin=varargin(2:end);
  elseif isequal(varargin{1},'-bar')
     lh=findall(gca,'type','hggroup'); % errorbar() handle
     bflag=1;
   % NB! errobar handles do have type hggroup (i.e. not of type line!)
  else lh=findall(gca,'type', 'line'); end

  if nargin>1 && isnumeric(varargin{1})
       cfac=varargin{1}; varargin=varargin(2:end);
  else cfac=[]; end

  getopt ('init', varargin);
     undo   = getopt ('undo'     ); if isempty(cfac)
     cfac   = getopt ('cfac',cfac); if isempty(cfac)
     cfac   = getopt ('fac', cfac); end, end
     cmin   = getopt ('cmin',  0.);
     back   = getopt ('2bg'      );
  varargin=getopt('get_remaining'); narg=length(varargin);

  if isempty(cfac), cfac=0.79; end

  if bflag
   % delete(findall(gca,'tag','blurl_bar'));
     for i=1:length(lh), ch=get(lh(i),'Children');
        c=get(lh(i),'Color'); set(ch,'Color',c,'ZData',[],'Visible','on');
        if undo, continue; end

        if cfac<1 && cfac>=0
           c = (1-cfac) * (c-1) + 1;
           for j=2 % 1:length(ch)
              x=get(ch(j),'XData');
              z=repmat(-0.01,size(x)); z(find(isnan(x)))=nan;
              set(ch(j),'ZData',z,'Color',c);
           end
        else
           set(ch(2),'Visible','off');
        end
     end

   % errorbars will still overlap with data points even though put to
   % background by setting z<0 above.
   % however, replotting the data without error bars still will find
   % itself partly covered with errorbars !-93_()*@!&$_
   % Wb,May07,09
   
   % if ~undo, hold on
   %    for i=1:length(lh), lo=lineopts(lh(i));
   %       x=get(lh(i),'XData');
   %       y=get(lh(i),'YData'); z=repmat(0.1,size(x));
   %       plot3(x,y,z,lo{:},'tag','blurl_bar');
   %    end
   % end

     if ~nargout, clear lh; end
     return
  end

  if undo
     for h=reshape(lh,1,[])
        u=get(h,'UserData');
        if isstruct(u) && isfield(u,'w') && isfield(u,'col')
        set(h,'LineWidth', u.w, 'Color', u.col); end
     end
  else
     cmin=sqrt(3)-cmin;
     for h=reshape(lh,1,[])
        w=get(h,'LineWidth'); c=get(h,'Color');
        set(h,'UserData',struct('w',w,'col',c));

        c2 = (1-cfac) * (c-1) + 1;
        if norm(c2)<cmin, c=c2; end
        set(h,'Color',c,'LineW',2);

        if ~isequal(h,'marker','none')
        set(h,'MarkerFaceColor',c,'MarkerEdgeColor',c); end
     end
  end

  if narg, set(lh,varargin{:}); end

  if back, mv2back(lh); end
  if nargout==0, clear lh; end

end

