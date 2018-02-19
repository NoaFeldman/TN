function gg=figdata(varargin)
% Function gg=figdata([opts])
%
%   Get data from existing figure which exists as an image
%   e.g. extracted from some pdf / publication.
%
% Options
%
%   '-x'  (re)specify x-range
%   '-y'  (re)specify y-range
%
% Wb,Nov30,08

% tags: pdf, png, xcoord, ycoord, points

  getopt('init',varargin);
     xflag=getopt('-x');
     yflag=getopt('-y');
     if getopt('-xy'), xflag=1; yflag=1; end
  varargin=getopt('get_remaining'); narg=length(varargin);

  if narg==1
    if ischar(varargin{1}) % assume image file name
        f=varargin{1}; varargin=varargin(2:end);
        if ~exist(f,'file'), wblog('ERR',...
        'no image file %s available',f); end

        fh=findall(groot,'type','figure','tag',mfilename);
        if ~isempty(fh), figure(fh); ah=getuser(fh,'ah'); else
           ah=smaxis(1,1,'tag',mfilename); setax(ah(1));
           A=imread(f); image(A);
        end

     elseif isscalar(varargin{1}) && isaxis(varargin{1})
        ah=varargin{1}; varargin=varargin(2:end);
        setax(ah);
     end
  end

  xsc=getuser(gcf,'xsc');
  if xflag || isempty(xsc) || ~isnumeric(xsc) || numel(xsc)~=3
     xsc=get_scale('x');
     setuser(gcf,'xsc',xsc);
  end

  ysc=getuser(gcf,'ysc');
  if yflag || isempty(ysc) || ~isnumeric(ysc) || numel(ysc)~=3
     ysc=get_scale('y'); 
     setuser(gcf,'ysc',ysc);
  end

  inl(1); k=0; gg={[]};

  while 1 % stop by pressing enter twice
     fprintf(1,['   collect xy-data (%d) by clicking to the figure '... 
     '(stop by pressing enter)...\n'],k);
     g=ginput; if isempty(g), inl(1); break; end

     k=k+1; gg{k}=g;
  end

  if length(gg)==1, gg=gg{1}; end

end

% -------------------------------------------------------------------- %
% -------------------------------------------------------------------- %

function sc=get_scale(x) 
% x = 'x' or 'y'

  while 1
     fprintf(1,...
     '   Determine %c-scale by clicking range [%c1 %c2] ...\n',x,x,x);
     figure(gcf); g=ginput;
     if size(g,1)==2, break; else
     fprintf(1,'   ERR exactly two points required!'); end
  end

  while 1
   % try
       i=input(sprintf(...
       '   Enter values [%c1 %c2] (including brackets): ',x,x));
   % catch, keyboard; continue; end

     if isnumeric(i) && numel(i)==2, break; else
     fprintf(1,'   ERR two values required!'); end
  end
  inl 1

  if x=='x'
     sc=[ g(1,1), (i(2)-i(1))/diff(g(:,1)), i(1) ];

     xl=get(gca,'XLim');
     xl= (xl-sc(1))*sc(2) + sc(3);

     for h=[ findall(gca,'type','line'); findall(gca,'type','image') ]'
        x=get(h,'XData'); x=(x-sc(1))*sc(2) + sc(3);
        set(h,'XData',x);
     end
     xlim(xl);
  else
     sc=[ g(1,2), (i(2)-i(1))/diff(g(:,2)), i(1) ];

     if ~isequal(x,'y'), error('Wb:ERR','invalid usage'); end
     yl=get(gca,'YLim');
     yl= (yl-sc(1))*sc(2) + sc(3);

     for h=[ findall(gca,'type','line'); findall(gca,'type','image') ]'
        y=get(h,'YData'); y=(y-sc(1))*sc(2) + sc(3);
        set(h,'YData',y);
     end
   % NB! images are drawn in axis ij mode
     set(gca,'YLim',sort(yl),'YDir','normal');
  end

end

% -------------------------------------------------------------------- %
% -------------------------------------------------------------------- %

