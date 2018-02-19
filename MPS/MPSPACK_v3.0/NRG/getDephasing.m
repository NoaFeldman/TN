function [it,iT,g0]=getDephasing(om,a0,Gamma,T,varargin)
% Function: [it,iT]=getDephasing(om,a0,Gamma,T [,OPTS])
% Options:
%
%   'a0fac'  factor applied to G to fulfill A(0)=pi·Gamma
%
%    Remainder of arguments (sucha as '-disc') are
%    handed over to the KramerKronig routine KKreal().
%
% See Micklitz et al, PRL 96, 226601 (2006)
% Wb,Sep05,07

% calc overall dephasing rate using -df/dw
% df/dw is symmetric w->-w; thus use |w| as this keeps numbers small
  ebo=exp(-abs(om)/T);
  dfw=(ebo/T) ./ (1+ebo).^2;

  getopt('init',varargin);
    a0fac =getopt('a0fac',1.);
    isdisc=getopt('-disc');
  varargin=getopt('get_remaining');

% NB! it involves different powers of G and subsequent integral
% involves sqrt(data) -> rather stick with continuous data!

  if isdisc
     a0=a0./diff2(om,'len');
  end

  if isreal(a0)
       g0=complex(KKreal(om,a0,varargin{:}),a0);
  else g0=a0; end

  if a0fac~=1, g0=g0/a0fac; end
  it=(pi*Gamma)*imag(g0) - (pi*Gamma)^2 * abs2(g0);

% integrated dephasing rate

% separate power ^p for real and imag part separately
% imag should not be there at the first place
% -> keep real real, and imag imag

% it1= it; it1(find(it1<0))=0; t1=it1; t1(find(t1==0))=Inf; t1=1./t1;
% it2=-it; it2(find(it2<0))=0; t2=it2; t2(find(t2==0))=Inf; t2=1./t2;

  dfw=dfw(:,ones(1,size(it,2)));

  iT=[
     pow_ri( intxy( om, dfw.*sqrt(it)),  2),... % d=3
   % pow_ri( intxy( om, dfw./sqrt(it)), -2)     % d=1 (numerically not very stable)
  ];

% if norm(iT)>1E2, keyboard, end

end

function zp=pow_ri(zz,p)

  iz=imag(zz); i=find(iz~=0); iz(i)=iz(i).^p;
  zp=complex( real(zz).^p, iz );
% zp=zz.^p;

end

% -------------------------------------------------------------------- %
% scratch
% -------------------------------------------------------------------- %

function scratch__()

% checking for negative 1/tau_phi (=it)
  ii=find(it<0);
  if ~isempty(it)
      eps=[ norm(it(ii))/norm(it), max(abs(it(ii)))/max(abs(it))];
      if all(eps<1E-6),     tag='NB!';
      elseif all(eps<1E-3), tag='WRN';
      else                  tag='ERR'; end

      wblog([tag ' neg. values for itau !?? (eps=%s)'],...
      mat2str2(eps,'fmt','%.3g','nofac'));

    % NB! sqrt(neg. values) strictly contributes to imag. part
    % -> keep for error measure upon analysis later
    % it(ii)=0;
  end

end

% -------------------------------------------------------------------- %
% -------------------------------------------------------------------- %

