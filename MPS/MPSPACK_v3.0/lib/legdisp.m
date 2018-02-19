function [lh,l2,ll]=legdisp(varargin)
% Function: [lh,l2,ll]=legdisp([opts])
%
%     display legend for line handles that carry an DisplayName.
%     returns handle to legend (lh), l2 returned line handles
%     by MatLab' slegend(), and ll actual line handles that appear
%     in this order in the legend.
%
% Options
%
%     '-flip'    flip order of entries
%     'pick',.   pick (sub)set in given order []
%     '-detach'  detach legend from axis (allows to create another legend)
%     '-erase'   white background but no border
%     '-eraset'  white background for text labels only
%     '-rm'      remove old legdisp axes sets
%     'xsc',..   scale legend horizontally (e.g. shrink if xsc<1)
%     'ysc',..   scale legend vertically (e.g. shrink if ysc<1)
%     'dx',...   shift of legend box / axis set
%     'dy',...   call to leg_mvtext if non-empty.
%     
%     'topt',..  further text options
%     'fs',..    font size of legend entries
%     
%     remaining options will be handed to as legend(...,OPTS).
%
%     NB! if -flip and pick are given,
%     then pick is done on the flipped handle set.
%
% Alternative usages
%
%  % get count and handles of line objects with non-empty disp only
%    [m,lh]=legdisp('-count');
%
% Wb,Oct05,05 ; Wb,Aug04,16

% by matlab>2013b legend is no longer of axis type,
% but is an object type of its own.
% for old version, see Archive/legdisp_160801.m
% => in case of detach, move to new axis handle

  if nargin==1 && isequal(varargin{1},'-count')
     l2=findall(gca,'type','line','-not','Disp','');
     lh=numel(l2); return
  end

  getopt('init',varargin);
   % usetag = getopt('-usetag');
     flip   = getopt('-flip'  ); % flip=~flip; % yet wrong order after repax(!)
     pick   = getopt('pick',{}); % correct order after plotting without flip!
     detach = getopt('-detach');
     rmold  = getopt('-rm');
   % lwidth = getopt('lwidth',[]);
     xsc    = getopt('xsc',[]);
     ysc    = getopt('ysc',[]);
     xl     = getopt('xl',[]);
     dx     = getopt('dx',[]);
     dy     = getopt('dy',[]);

     topt   = getopt('topt',{});
     fs     = getopt('fs',[]);

     if getopt('-erase'); erase=1;
     elseif getopt('-eraset'); erase=2; else erase=0; end

  varargin=getopt('get_remaining');

  if isempty(pick) && ~iscell(pick)
     if ~nargout, clear lh; else lh=[]; l2=[]; end
     return
  end

  ah=gca; o={'type','line','visible','on'};

  lh=flipud(findall(ah,o{:}));

% arrowsafe do not contain property field `DisplayName' (!?)
  h=findall(ah,o{:},'tag','arrowsafe');
  if ~isempty(h)
   % NB! setdiff() also sorts data in lh => changes order of handles!
   % Wb,Mar05,11
     [q,is]=setdiff(lh,h);
     lh=lh(sort(is)); % preserve original handle order!
  end

% arrowsafe do not contain property field `DisplayName' (!?)
  nn=get(lh,'DisplayName'); if length(lh)==1, nn={nn}; end
  n=length(lh); m=zeros(1,n);

  for i=1:n
     if ~isempty(nn{i}), m(i)=1;
      % when plotting multiple columns and using 'Disp','XXX'
      % DisplayName' is actually set to getcolumn(XXX,#)
      % -> replace to XXX(:,#)
        nn{i}=regexprep(nn{i},'^getcolumn\((.*),([0-9])\)','$1(:,$2)');
     end
  end

  ii=find(m); n=numel(ii);

  if flip, ii=fliplr(ii); end
  if ~isempty(pick)
     if ischar(pick), eval(sprintf('ii=ii(%s);',pick));
     else
        i=find(pick>n);
        if ~isempty(i), wblog(2,'WRN',...
          'legdisp() pick out of range (%g/%g; skip)',max(pick),n);
           pick(i)=[];
        end
        ii=ii(pick);
     end
     n=length(ii);
  end

  if ~n, legend off
     if ~nargout, clear lh; else lh=[]; l2=[]; end
     return
  end

  ll=lh(ii); nn=nn(ii);

% [lh,l2]=legend(...) appears to call drawnow in matlab/2016a
  [lh]=legend(ll,nn{:},varargin{:});

  if erase==1
   % tiny offset (otherwise frame still appears in black when printing)
     if isnumeric(lh) % matlab/2013b
          set(lh,'Color','w','XColor',[1 1 1]-1E-12,'YColor',[1 1 1]-1E-12)
     else set(lh,'Color','w','EdgeColor',[1 1 1]-1E-12);
     end
  else
     legend('boxoff');
  end

  if ~isempty(dx), mvaxis(lh,dx); end

% remove existing legdisp sets that share exactly the same position
  if rmold
     h0=findall(groot,'Type','axes','Position',get(lh,'Position'));
     if numel(h0)
        if 1, delete(h0);
           t=sprintf(',%.2f',get(lh,'Position'));
           t=sprintf('legdisp:lh[%s]',t(2:end));
           wblog(' * ','got %g legend(s) at %s!',numel(h0),t);
        end
     end
  end

  if ~detach && erase<=1 && ...
     isempty(xl) && isempty(dy) && isempty(xsc) && isempty(ysc)
     if ~isempty(fs), set(lh,'FontSize',fs); end
     setax(ah); if ~nargout, clear lh; end
     return
  end

% -------------------------------------------------------------------- %
% switch legend to axes handle // DETACH
% since for matlab>2013a legend is a type of its own // Wb,Aug04,16
% -------------------------------------------------------------------- %

  lh_=lh; th=zeros(1,n);
  lh=axes('Position',get(lh_,'Position'));

  if isempty(fs), fs=get(lh_,'FontSize'); end
  topt={topt{:},'FontSize',fs,'Tag','leg:text'};

% redraw content of legend (it appeares difficult to extract
% line+text data from legend for matlab>2013b)

  los={
     'LineStyle','LineWidth','Color',...
     'Marker','MarkerSize','MarkerEdgeColor','MarkerFaceColor'
  };

  to={'Tag','leg:line'};

  for i=1:n, h=ll(i); lox=get(h,los); y=n-(i-1);
     q={los{1:3}; lox{1:3}}; plot([-2 0],[y y],q{:},to{:});
     if i==1, hold on
        axis([-2.5 5, 0.5 n+0.5]); % first guess (see update below)
     end
     if ~isequal(lox{4},'none')
        q={los{3:end}; lox{3:end}}; plot(-1,y,q{:},to{:});
     end
   % use tex-format for all entries (otherwise there is slight
   % difference in font between tex-entries and plain text
     th(i)=text(0.5,y,[nn{i} '^{}'],topt{:});
  end

  xy=get(th,'Extent'); if iscell(xy), xy=cat(1,xy{:}); end
  axis([-2.5, 0.5+1.2*max(xy(:,3)), [1 n]+[-.5 .5]*max(xy(:,4))]);

% link to corresponding axis set
  setuser(lh,'Parent',ah);

  if 0
     t=sprintf(',%.2f',get(ah,'Position'));
     t=sprintf('legdisp:ah[%s]',t(2:end));

     h0=findall(groot,'Type','axes','Tag',t);
     if ~isempty(h0), t
        wblog('WRN','got %g user legend(s) for given axis set',numel(h0));
     end
     set(lh,'Tag',t);
  end

% tiny offset (otherwise frame still appears in black when printing)
  if erase==1
     setax(lh); box on
     set(lh,'XColor',[1 1 1]-1E-12,'YColor',[1 1 1]-1E-12,'Visible','on',...
     'XTick',[],'YTick',[])
  else
     if erase, set(th,'BackgroundC','w','Margin',1E-3); end
     axis off
  end

% setax(ah); return % TST Wb,Aug04,16 

% adjust y-position of legend entries (applied to text)
  if ~isempty(dy), leg_mvtext(lh,dy); else
     for i=1:n, dy=getuser(ll(i),'leg_dy');
      % alternatively, this may already be set in caller routine
        if ~isempty(dy), h=th(i);
           if numel(dy)==1, dy=[0 dy 0]; else dy(end+1:3)=0; end
           set(h,'Pos',get(h,'Pos')+dy);
        end
     end
  end

% NB! requires -detach (!) // Wb,Dec04,10
% otherwise, upon redrawing axis set, legend becomes strongly distorted!
  if ~isempty(xsc) && ~isequal(xsc,1)
     if numel(xsc)==2, set(th,'Units','data');
        if iscell(xsc)
             set(lh,'XLim',xsc{2}); xsc=xsc{1};
        else set(lh,'XLim',[0 1/xsc(2)]); end
        set(th,'Unit','norm')
     end
   % x=get(lh,'XLim'); set(lh,'XLim',mean(x)+diff(x)*[-0.5 0.5]/xsc(1));
     p=get(lh,'Pos'); x=(p(1)+p(3)/2) + [-0.5 0.5]*(p(3)*xsc(1));
     set(lh,'Pos',[x(1), p(2), diff(x), p(4)]);
   % set(lh,'XColor','r','YColor','r');
   % get(lh,'XLim')
  end

% NB! requires -detach (!)
% otherwise, upon redrawing axis set, legend becomes strongly distorted!
  if ~isempty(ysc) && ~isequal(ysc,1)
     if numel(ysc)==2, set(lh,'YLim',0.5+[-1 1]*(0.5/ysc(2))); end
     set(th,'Units','data');
     p=get(lh,'Pos'); y=(p(2)+p(4)/2) + [-0.5 0.5]*(p(4)*ysc);
     set(lh,'Pos',[p(1), y(1), p(3), diff(y)]);
     set(th,'Unit','norm')
  end

  if ~isempty(xl), set(lh,'Xlim',xl); end

% NB! requires -detach (!?)
% if ~isempty(lwidth), set(findall(lh,'type','line'),'LineW',lwidth); end

  setax(ah); legend hide
  if ~nargout, clear lh; end

end

