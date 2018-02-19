function [h,th]=ymark(varargin)
% usage: ymark(ypos [,opts])
% mark y position in current figure
%
% Options:
%
%    'istr',str  info string positioned by default at relative xpos=0.12
%    'istr',{str, [xpos,][text options]}
%
%    '-rm'   remove all existing ymarks
%    '-fix'  redraw ymarks so they stretch over full xlim-range
%    '-bg'   move ymark into background
%    'xl',.. limits of horizontal marker in normalized units
%
%    remainder of options are used as options to ymark line to be drawn.
%
% See also xmark, ypatch.

  getopt('init',varargin);
     istr =getopt('istr','');
     rmall=getopt('-rm'  );   % remove all xmarks
     fxall=getopt('-fix');    % fix length of markers assuming ylim changed
     bflag=getopt('-bg');     % move to background
     fgflg=getopt('-fg');     % keep in foreground
     xl1  =getopt('xl',[]);
  varargin=getopt('get_remaining'); narg=numel(varargin);

  if rmall
     delete(findall(gca,'tag','ymark'));
     if ~narg, return; end
  end

  xl=xlim; gotlog=isequal(get(gca,'XScale'),'log');
  if gotlog && xl(1)==0
     t=get(gca,'XTick'); if length(t>1), xl(1)=t(1);
     else wblog('WRN','xlim returns 0 as left boundary'); end
  end

  if ~isempty(xl1)
     if numel(xl1)~=2, error('Wb:ERR','invalid xl(im) values'); end
     if gotlog, xl=log(xl); end
     xl=xl(1)+xl1*diff(xl);
     if gotlog, xl=exp(xl); end
  end

  if fxall
     for m=findall(gca,'type','line','tag','ymark')', set(m,'XData',xl); end
     if ~narg, return; end
  end

  if length(varargin)==0, return; end
  if isnumeric(varargin{1})
     y=varargin{1}; varargin=varargin(2:end);
  end

  hold on

  if mod(length(varargin),2)
       fmt=varargin(1); varargin=varargin(2:end);
  else fmt={'Color', [1  1 .7],'LineW',4}; bflag=1; end

  if fgflg, bflag=0; end

  h=plot(xl, [y(:) y(:)]', fmt{:},'tag','ymark');
  if bflag, mv2back(h); end

  if ~isempty(varargin)
  set(h,varargin{:}); end

  if ~isempty(istr)
     xt=0.12; n=numel(y); th=zeros(n,1);

     if iscell(istr), topt=istr(2:end); istr=istr{1};
        if length(topt) && isnumeric(topt{1})
        xt=topt{1}; topt=topt(2:end); end
     else topt={}; end

     c=get(h(1),'Color');
     if norm(c)>0.8, c=0.8*c/norm(c); end

     for i=1:n
      % no xt here (bad for log-plots)
        th(i)=text(1,y(i),istr,...
        'Color',c,'FontWeight','bold','tag','ymark',topt{:});
      % NB! y-axis may be log-plot(!)
        set(th(i),'Units','norm'); p=get(th(i),'Pos'); p(1)=xt;
        set(th(i),'Pos',p,'Units','data');
     end
  end

  if ~nargout, clear h; end

end

