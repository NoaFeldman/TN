function th=show_deg(Ek,varargin)
% Function th=show_deg(Ek [,ah,opts])
%
%    show degeneracies of line plot (e.g. NRG flow diagram).
%
% Options
%
%    'dz',..    internal multiplet dimension (one row for each value Ek)
%    'qq',..    corresponding symmetry labels (one row for each value Ek)
%
%    'fs',..    fontsize of labels (8)
%    'tag',...  tag to apply to text labels (':deg:')
%    'fs',..    fontsize of labels (8)
%    'eps',...  consider levels degenerate within eps (1E-3)
%    'x0',...   x-position of labels (if 0<x0<1, assume relative x-axis position)
%    'xt',...   same as x0, but always in absolute units
%    'dx',...   x-offset if non-degenerate y-values by less than dy
%    'dy',...   threshold to use dx (see dx)
%    'n',...    maximum number of values to consider in sorted Ek
%    'y2',...   ignore values Ek>y2 (default: 0.6 of ylim range)
%   
%    'nmax', .. do not show if number of closeby levels exceeds nmax (Inf)
%    'lmax', .. show at most the lowest lmax levels for each symmetry (Inf)
%
% see also nrg_plot.m
% Wb,Mar18,08

  if nargin>1 && isaxis(varargin{1})
     setax(varargin{1}); varargin=varargin(2:end);
  end

  fs=8;
  getopt('init',varargin);

     dz   = getopt('dz',[]); % Wb,Apr14,12
     qq   = getopt('qq',[]); % if ~isempty(dz) || ~isempty(qq), fs=6; end
     vflag= getopt('-v');

     tag  = getopt('tag',':deg:');
     fs   = getopt('fs',fs);
     xt   = getopt('xt',[]); if isempty(xt)
     x0   = getopt('x0',0.8); end
     dx   = getopt('dx',[]);
     dy   = getopt('dy',0.25);
     nmax = getopt('nmax',Inf);
     lmax = getopt('lmax',Inf);
     y2   = getopt('y2',[]);
     eps  = getopt('eps',1E-3);
     n    = getopt('n',128);
  varargin= getopt('get_remaining'); narg=length(varargin);

  if isempty(qq)
       va={'VerticalAl','bottom','Margin',1E-3}; % ,'BackgroundC','w'
  else va={'VerticalAl','bottom'}; end
  topt={'HorizontalAl','center',va{:},'tag',tag,'FontSize',fs,varargin{:}};
  % 'VerticalAl','middle'

  yl=ylim; if isempty(y2), y2=yl(1)+0.6*diff(yl); end
  n=numel(Ek);

  if ~isvector(Ek)
     error('Wb:ERR','\n   ERR invalid usage (expecting vector for Ek)'); end
  if ~isempty(dz) && size(dz,1)~=n, whos dz Ek
     error('Wb:ERR','\n   ERR invalid usage (size mismatch of dz)'); end
  if ~isempty(qq) && size(qq,1)~=n, whos qq Ek
     error('Wb:ERR','\n   ERR invalid usage (size mismatch of qq)'); end

  [Ek,is]=sort(reshape(Ek,[],1));
  Ek(find(isnan(Ek)))=Inf;

  if ~isempty(dz), dz=dz(is,:); end
  if ~isempty(qq), qq=qq(is,:); end

% i1=1; I=i1; n=min(n,length(Ek));
% while i1<n && Ek(i1)<y2
%    i1=i1+find(Ek(i1+1:n)>(Ek(i1)+3*eps),1); if isempty(i1), break; end
%    I(end+1)=i1;
% end
  I=[1; 1+find(diff(Ek)>=eps) ]; % , n+1
  i=find(Ek(I)>y2); I(i(2:end))=[];

  xl=xlim; % disp([xl x0])
  if ~isempty(xt), x0=xt;
  elseif x0<=1, x0=xl(1)+x0*diff(xl); end

  if isempty(dx), dx=diff(xl)/50; end
  m=length(I); xi=repmat(x0,m,1);

% if x0==0, keyboard, end

  Ex=Ek(I); i1=1; n2=numel(Ex); nx=numel(dx); count=0;
  while i1<m
     i2=i1+find(Ex(i1+1:n2)>(Ex(i1)+dy),1); if isempty(i2), break; end
     count=count+1; if count>lmax, xi(i1:end)=Inf; break; end

     if i2>i1+1, q=(i2-i1-1)/2; q=(-q:q);
        if nx==1, q=q*dx;
        else
         % expect range dx = [ dx_min, dx_max ] // Wb,Apr18,11
         % got set of 2: use dx_max
         % got set of 3: use half-way between dx_min and dx_max
         % got set of more than 3: use dx_min
           nq=numel(q); if nq>nmax, q(:)=Inf;
           else
              if nq==2, q=q*dx(2);
              elseif nq==3, q=q*mean(dx(1:2));
              else q=q*dx(1); end
           end
        end
        xi(i1:(i2-1))=x0+q;
     end
     i1=i2;
  end

% if x0==0, keyboard, end

  de=0.1*mean(diff(Ex));

  for k=2:m, i=I(k-1):(I(k)-1); d=numel(i);
     if isempty(dz)
        dstr=sprintf('%g',d);
     else
        g=prod(dz(i,:),2);
        if numel(g)>1 && vflag % && isempty(qq)
           s=sprintf('%+g',g);
           dstr=sprintf('%s=%g',s(2:end),sum(g));
        else dstr=sprintf('%g',sum(g)); end
     end
     th(k-1,1)=text(xi(k-1),mean(Ek(i)), dstr,topt{:},'VerticalA','middle');

     if ~isempty(qq)
        if isempty(dz)
           th(k-1,2)=text(xi(k-1),mean(Ek(i)), ...
              mat2str2(qq(i,:),'fmt','%2g'), ...
              topt{:},'VerticalA','top','FontS',fs-2 ...
           );
        else
           fmt=[ '[ ' repmat('%2g ',1,size(qq,2)) '] (%g)\n' ]; g=prod(dz(i,:),2);
           th(k-1,2)=text(xi(k-1),mean(Ek(i)), ...
              sprintf(fmt,[qq(i,:),g]'),...
              topt{:},'VerticalA','top','FontS',fs-2 ...
           );
        end
     end
  end

  if ~nargout, clear th; end

end

