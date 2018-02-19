function [ff,N,I] = getNRGcoupling(Gamma,Lambda,N,varargin)
% Funciton: [ff,N,I] = getNRGcoupling(Gamma,Lambda,N [,OPTS])
%
%   Gamma  coupling of bath states to impurity
%   Lambda NRG descretization parameter
%   N      total length of Wilson chain, i.e. length(ff)=N-1
%
% Options
%
%  'z',..  Oliveira shift in log. disrcretization (0)
%  'r',..  pseudo-gap model with rho~omega^r (0)
%
%  '-x'    extend to energy scales below double precision by fitting
%  '-w'    using standard Wilson discretization (i.e. not -AL)
%  '-AL'   (deprecated) accounting for wavefunction renormalization
%          A_Lambda originating from discretization using corrected
%          energies in line with Oliveria (2005), Rok Zitko (2009).
%
% Wb,May23,06

  if nargin<3
     eval(['help ' mfilename]);
     error('Wb:ERR','Invalid number of input arguments');
  end

  if nargout>1, I=[]; end

% for compatibility reasons
  if nargin>3 && isnumeric(varargin{1})
     z=varargin{1};
     varargin=varargin(2:end);
  else z=[]; end

  oli={'-AL'};

  getopt ('init', varargin);
    z = getopt('z',z);
    r = getopt('r',0);
    xflag = getopt('-x');

    if getopt('-p'); qflag=-1; else qflag=getopt('-q'); end

    if getopt('-w'), oli={};
    elseif getopt('-AL')
       oli={'-AL'};
       wblog('WRN','-AL is deprecated'); % since -AL is default setting!
    end
  varargin=getopt('get_remaining'); oli={ oli{:}, varargin{:} };

  if isempty(z)
     global OLIZ % Oliveira Z
     if isempty(OLIZ), z=0; else z=OLIZ; end
  end

  if N<=1, ff=[]; return; end

  if ~isnumber(Gamma,Lambda,N)
     error('Wb:ERR','invalid usage (need scalar arguments)'); end

  if xflag
     n=ceil(-2*log(1E-15)/log(Lambda));
     setopts(oli,'-q',{'Nmax',-1});
     if N<n+10, xflag=0;
     else N_=N; N=n; end
  end

  if ~isempty(r) && r~=0
     ff=oliveira_vdisc(Gamma,Lambda,'z',z,'r',r,'N',N,'nolog');
     return
  end

  Delta=1/Lambda;

  if ~qflag, oli{end+1}='-q';
  elseif qflag>1, oli{end+1}='-Q'; end

% [alpha,beta]=oliveira(Gamma,Lambda,z,'nolog');
  [fx,ex,I]=oliveira(Gamma,Lambda,z,'N',N,oli{:});
     I.fx=fx;
     I.ex=ex;

  if xflag
     wblog(' * ','extrapolating NRG couplings (%d->%d)',N,N_);
     [fx,N]=fix_couplings(fx,N_-1);
  elseif N-1<=numel(fx), fx=fx(1:N-1); ex=ex(1:N-1);
  else
     n=N; N=numel(fx)+1; wblog('WRN',...
    'Wilson chain will be truncated N=%g->%g (%.3g)',n,N,fx(end));
  end

  if z==0 && isempty(oli)
     ii=0:N-3;

     xi = (1-Delta.^(ii+1))./sqrt((1-Delta.^(2*ii+1)).*(1-Delta.^(2*ii+3)));
   % in the limit of large i, with Delta<1 => xi=1.

     ff = [ sqrt(2*Gamma/pi), (1+Delta)/2*(Delta.^(ii/2)) .* xi ];
   % ff = [ sqrt(2*Gamma/pi), (1+Delta)/2                 .* xi ];

     I.fx=ff;
     if max(abs(ff-fx))>1E-12
     error('Wb:ERR','Oliveira - mismatch in couplings.'); end
  else
     ff=fx;
  end

end

% try to extend to energy scales below double precision by fitting
% while respecting even/odd behavior % Wb,Aug31,12

function [fx,N]=fix_couplings(fx,N);

  n=length(fx); i0=ceil(0.66*n);

  vflag=0;
  if vflag>1
     wblog('TST','applying even/odd contrast (factor 20)'); 
     fx(1:2:end)=20*fx(1:2:end);
  end

  for k=1:2
     i=i0+k:2:n; y=log(fx(i));
     p{k}=polyfit(i,y,2); e(k)=norm(y-polyval(p{k},i));
  end

% NB! e<<1 corresponds to relative error
% 1-y/y0 = |1-exp(log(y)-log(y0))| \simeq |log(y)-log(y0)|
  if norm(e)<1E-4
     for k=1:2
        i=i0+k:2:N; i=i(find(i>n));
        fx(i)=exp(polyval(p{k},i));
     end
     wblog(iff(norm(e)<1E-8,' * ','WRN'), [
      'extending couplings L=%d (%.3g) => L=%d (%.3g)\n' ...
      '(at rel. accuracy %.3g)'], n,fx(n),N,fx(N),norm(e)); 
  end

  if ~vflag, return, end

ah=smaxis(1,1,'tag',mfilename); header('%M'); addt2fig Wb
setax(ah(1,1))

  semilogy(fx,'o-'); sms(4); hold on
  i=n+1:N; h=semilogy(i,fx(i),'ro'); sms(h,4);

  xmark([i0,n,N],'k:'); xlim([0 N+1]);

% keyboard

end

