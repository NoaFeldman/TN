function [h,th]=xmark(varargin)
% usage: xmark(xpos, string [,color])
% mark x position in fig.
%
% Options
%
%   'istr',..   display istr with last xmark line
%   'itop'      display istr on upper end of xmark line (default: bottom)
%   'iTOP'      display istr on top (outside axis set)
%   'top'       keep xmark lines front (default: moved to background)
%   'yl',[]     ylim of xmark lines in normalized units
%   '-rm'       remove all marks
%
% Example: put 'istr' at ypos=0.8 (normalized) in vertical orientation
%
%   'istr',{'t_{Kondo}',0.8,'Rotation',90})
%
% See also xpatch.
% Wb,Mar07

  getopt('INIT',varargin);
     istr =getopt('istr','');
     itop =getopt('itop',''); if ~isempty(itop), istr=itop; itop=1; end
     iTOP =getopt('iTOP',''); if ~isempty(iTOP), istr=iTOP; itop=2; end
     top  =getopt({'-fg','top'});
     yl   =getopt('yl',[]); % ypos in normalized units
     rmall=getopt('-rm');   % remove all xmarks
     gflag=getopt('-g' );   % gray solid lines
     fxall=getopt('-fix');  % fix length of markers assuming ylim changed
     r90  =getopt('-90' );
  varargin=getopt('get_remaining');

  if rmall
     delete(findall(gca,'tag','xmark')); % 'type','line'
     if nargin<2, return; end
  end

  dm=0.005;

  if fxall
     yl=ylim; s=get(gca,'YScale');
     if isequal(s,'log'), if yl(1)==0, yl(1)=1E-16; end % safeguard
          dy=(yl(2)/yl(1))^dm; yl=[yl(1)*dy, yl(2)/dy];
     else dy=dm*diff(yl); yl=[yl(1)+dy, yl(2)-dy]; end

     hh=findall(gca,'type','line','tag','xmark')';
     for m=hh, set(m,'YData',yl); end
   % if ~isempty(hh), wblog('<i>','fixed %g xmarks.',length(hh)); end
     if nargin<2, return; end
  end

  if length(varargin)==0, return; end
  if isnumeric(varargin{1})
     x=varargin{1};
     varargin=varargin(2:end);
  else x=[]; end

  hold on

  if mod(length(varargin),2)
       fmt=varargin(1); varargin=varargin(2:end);
  else fmt={'Color', [1  1 .7],'LineW',4}; end

  if isempty(yl)
     yl=ylim; s=get(gca,'YScale');
     if isequal(s,'log')
          dy=(yl(2)/yl(1))^dm; yl=[yl(1)*dy, yl(2)/dy];
     else dy=dm*diff(yl); yl=[yl(1)+dy, yl(2)-dy]; end
  else l=ylim;
     if isequal(get(gca,'YScale'),'log') && all(l>0)
          yl=l(1)*(l(2)/l(1)).^yl;
     else yl=l(1)+diff(l).*yl; end
  end

  if isempty(x), return; end

  if any(isnan(yl))
   % wblog('WRN','xmark - invalid ylim [(%.4g %.4g]',yl);
     yl=getylim('-data');
  end

  h=plot([x(:) x(:)]', yl,fmt{:},'tag','xmark');
  if ~top, mv2back(h); end

  if gflag
  set(h,'Color',[1 1 1]*.9, 'LineW',0.5); end

  if ~isempty(varargin)
  set(h,varargin{:}); end

  if ~isempty(istr)
     c=get(h(1),'Color'); fw={'FontWeight','bold'};
     if norm(c)>0.8, c=0.8*c/norm(c); end

     ha='center'; yt=0.12;
     if r90, rot={'Rotation',90}; else rot={}; end

     if iscell(istr), topt=istr(2:end); istr=istr{1};
        if length(topt) && isnumeric(topt{1})
        yt=topt{1}; topt=topt(2:end); end
     else topt={}; end

     switch itop
        case 1, yt=0.90; 
        case 2, yt=1.15; c='k'; fw={};
     end

     n=numel(x); th=zeros(n,1);

     xl=xlim; if x(end)>=xl(1) && x(end)<=xl(2)
      % no yt here (bad for log-plots)
        th=text(x(end),1,istr, rot{:},'Color',c,fw{:},...
          'tag','xmark','UserD',yt,...
          'HorizontalAl',ha,'BackgroundC','w',topt{:});
      % NB! y-axis may be log-plot(!)
        set(th,'Units','norm'); p=get(th,'Pos'); p(2)=yt;
        set(th,'Pos',p); set(th,'Units','data');
     end
  end

  if ~nargout, clear h; end

end

