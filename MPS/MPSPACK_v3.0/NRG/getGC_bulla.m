function [ac,I]=getGC_bulla(varargin)
% Function: [ac,I]=getGC_bulla(om, gi,fi [,opts])
% Function: [ac,I]=getGC_bulla([om gi], [om fi] [,opts])
%
%    calculate improved spectral function using Bulla'98
%    with gi and fi preferentially already but not necessarily
%    smoothened (note that KramersKronig within discrete data
%    is also possible, but when smoothening this deteriorates
%    somewhat the sum-rule for spectral function).
%
%    (SIAM) parameters are read from global/base `param'.
%
% Options
%
%   'omfac',..   linear interpolation for energies < min(om)*omfac (2)
%
%   NB! all remaining options are forwarded to KKreal.
%
% Wb,Oct22,05 ; Wb,Jan23,11

% outsourced from $VER//xGF

  for narg=1:nargin
     if ~isnumeric(varargin{narg}), narg=narg-1; break; end
  end

  if ~nargin || narg<2 || narg>3
     fprintf(1,'\n'); eval(['help ' mfilename]);
     if nargin || nargout, error('Wb:ERR','invalid usage'), end, return
  end

  if narg==2
   % getGC_bulla([om,gi],[om,fi] [,opts])
     [gi,fi] = deal(varargin{1:2}); varargin=varargin(3:end);
      % with gi = [om,gg], fi = [om,ff]
      % => 1st column in gi and fi is omega data (checked to be the same)
      % => remaining columns are actual gg / ff data
     om=gi(:,1);
        if length(om)~=size(fi,1) ||  norm(om-fi(:,1))>1E-15
        error('Wb:ERR','ERR gi and fi must have the same omega range!'); end
     gi=gi(:,2:end); fi=fi(:,2:end);

  elseif narg==3,
   % getGC_bulla(om,gi,fi [,opts])
     [om,gi,fi] = deal(varargin{1:3}); varargin=varargin(4:end);

     if ~isvector(om) || length(om)~=size(gi,1) || length(om)~=size(fi,1)
     error('Wb:ERR','invalid omega (size mismatch)'); end

  else error('Wb:ERR','invalid usage'); end

  getopt('init',varargin);
     omfac=getopt('omfac',2);
     iBflag=getopt('-iB');
  varargin=getopt('get_remaining');

  if ~isvector(om)
     eval(['help ' mfilename]);
     if nargin || nargout, error('Wb:ERR','invalid usage'), end, return
  else om=om(:); end

  global param
  if isempty(param), getbase param; end
  if isempty(param), error('Wb:ERR','ERR failed to get param(eters)'); end

  if isfield(param,'Gamma') && isfield(param,'epsd') && isfield(param,'U')
     B=0; if isfield(param,'B') B=param.B; end
     if iBflag, B=-B; end

     [gc,US,I]=getBulla_SIAM(om,gi,fi,...
     param.U,param.epsd,param.Gamma,B,iBflag,varargin);

  elseif isfield(param,'Gamma') && isfield(param,'JH') && ...
    ~isfield(param,'U') && ~isfield(param,'epsd')

   % B=0; if isfield(param,'B') B=param.B; end
   % [gc,US,I]=getBulla_JH(om,gi,fi, param.JH,Gamma,B, varargin);
     error('Wb:ERR','use getGC_JJH.m instead');

  else
     error('Wb:ERR','ERR do''nt know what to do with param(eter) set');
  end

  add2struct(I,gc,US);

  if omfac>0
   % linear interpolation for range around om==0, Wb,Jun16,09

     aom=abs(om); x0=min(aom)*omfac;

        i0=find(aom<x0);
        i1=find(om==max(om(find(om>=-x0))),1);
        i2=find(om==min(om(find(om>=+x0))),1);

        add2struct(I,omfac,x0);

     if isempty(i0) || isempty(i1) || isempty(i2), return; end
     om0=om(i0); x12=[om(i1),om(i2)];

     for k=1:size(fi,2)
      % NB! polyfit/val accept complex y-data
        gc(i0,k) = polyval(polyfit(x12,[gc(i1,k),gc(i2,k)],1),om0);
        US(i0,k) = polyval(polyfit(x12,[US(i1,k),US(i2,k)],1),om0);
     end
  end

  ac=-(1/pi)*imag(gc);

end

% -------------------------------------------------------------------- %
% -------------------------------------------------------------------- %

function [gc,US,I]=getBulla_SIAM(om,gi,fi, U,epsd,Gamma,B,iBflag, oKK)

  m=size(fi,2);
  if m>2, error('Wb:ERR','invalid input data (%g)',size(fi,2)); end
  if B~=0 && m<2 && ~iBflag
     wblog('WRN','assuming spin up for improved spectral function');
  end

  e0=epsd+[-B/2, +B/2]; % e0=e0(1:m);

  I=add2struct('-',U,epsd,Gamma,B,e0);

  dd = get_delta(om,Gamma);

% NB! include factor A(w) = -(1/pi) Im G(w) => Im G(w) = -pi * A(w)
% same factor of ff => Uf actually independent of factor applied to gi AND fi
  gi=(-pi)*gi;
  fi=(-pi)*fi;

  gg=zeros(size(gi));
  ff=zeros(size(fi));

  for k=1:m
     gg(:,k) = complex( KKreal(om,gi(:,k),oKK{:}), gi(:,k) );
     ff(:,k) = complex( KKreal(om,fi(:,k),oKK{:}), fi(:,k) );

     Uf = U*(ff(:,k)./gg(:,k));

     gc(:,k) = 1./( om - e0(k) - dd - Uf);
     US(:,k) = Uf;
  end

  add2struct(I,gg,ff);

% keyboard

end

% -------------------------------------------------------------------- %
% -------------------------------------------------------------------- %

function [gc,gs,I]=getBulla_JH__(om,gi,fi, JH,Gamma,B, oKK)
% deprecated; see getGC_JJH.m

  m=size(fi,2);
  if m>2, error('Wb:ERR','invalid input data (%g)',size(fi,2)); end
  if B~=0 && m~=2, wblog('WRN',...
  'assuming spin up for improved spectral function'); end

  e0=[-B/2, +B/2]; I=add2struct('-',JH,Gamma,B,e0);

  dd = get_delta(om,Gamma);

  for k=1:m
     gg = complex( KKreal(om,gi(:,k),oKK{:}), gi(:,k) );
     ff = complex( KKreal(om,fi(:,k),oKK{:}), fi(:,k) );

     ss(:,k) = ff./gg;
     gc(:,k) = -1./( om - e0(k) - dd - JH*ss(:,k));
  end

end

% -------------------------------------------------------------------- %
% -------------------------------------------------------------------- %

function dd=get_delta(oz, Gamma)

  if isreal(oz), oz=complex(oz, +1E-12); end % causal Greens function

  dd = (Gamma/pi)*log((oz+1)./(oz-1));

end

% -------------------------------------------------------------------- %
% -------------------------------------------------------------------- %

