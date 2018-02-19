function c=getcolor(varargin)
% Function: getcolor(i,[CC])
%
%    get i-th color of color order of current axis set
%    taken modulo if i exceeds colormap range
%
% Options
%
%    CC  user defined color order
%
% See also recolor.m Wb,Apr20,07
% Wb,Jan18,08

  if nargin==1, q=varargin{1};
     if numel(q)==3 && all(q>=0) && all(q<=1), c=q; return; end
     if ischar(q) && numel(q)==1
        if ~isempty(find('brmcyk'==q)), c=q; return; end
        if isequal(q,'g'), c=[0 0.5 0]; return; end
     end
  end

  getopt('init',varargin);
      flip=getopt('-flip');
      cmap=getopt('-cm'  ); % used colormap instead of colororder
  varargin=getopt('get_remaining');
  narg=length(varargin);

  if narg && ischar(varargin{1})
       uflag=1; ucol=varargin{1}; varargin=varargin(2:end);
  else uflag=0; end

  narg=numel(varargin);
  if uflag && narg || ~uflag && (narg<1 || narg>2)
     eval(['help ' mfilename]);
     if nargin || nargout, error('Wb:ERR','invalid usage'); end
     return
  end

  if uflag
     if ~isempty( regexp(ucol,'^gray[0-9][0-9]$') )
        c=repmat(0.01*str2num(ucol(5:end)),1,3);
     else
        switch ucol
           case 'brown',     c = [ 0.88  0.68  0.3  ];

           case 'blue',      c = [ 0.1   0.5   1    ];
           case 'lightblue', c = [ 0.5   0.75  1    ];

           case 'orange',    c = [ 1     0.5   0    ];
           case 'orange.in', c = [ 0.88  0.68  0.34 ];
         % case 'gray80',    c = [ 0.80  0.80  0.80 ]; see above

         % accept standard color specs (eg. see usage of *this in getcolor())
           case 'g',         c = [ 0 .5 0 ];
           case 'y2',        c = [ .85 .85 0 ];

        otherwise
           disp(ucol)
           error('Wb:ERR','unknown user color specification');
        end
     end
     return
  end

% alternative usage: index [,colormap]

  i=varargin{1};
  if isequal(size(i),[1 3]) && all(i<=1) && all(i>=0)
   % accept standard color specs (eg. see usage of *this in getcolor())
     c=i; return
  end

  if narg<2
     if cmap, CC=colormap;
     else CC=get(gca,'ColorOrder'); end
  else
     CC=varargin{2};
     if ~isnumeric(CC) || size(CC,2)~=3
     error('Wb:ERR','invalid color set'); end
  end

  nc=size(CC,1); i=i(:);

  if flip
       j=mod(nc-i,nc)+1;
  else j=mod(i-1, nc)+1; end

  if isequal(j,round(j))
     c=CC(j,:);
  else
     j1=floor(j); j2=ceil(j); j2(find(j2>nc))=nc;
     x=repmat(j2-j,1,3);

     c=CC(j1,:).*x + CC(j2,:).*(1-x);
  end

end

