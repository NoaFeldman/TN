function [b,a,I] = oliveira(Gamma, Lambda, z, varargin)
% Function: oliveira
% Logarithmic discretization using Oliveira's trick (PRB '90)
%
% Usage: [beta, alpha] = oliveira(Gamma, Lambda, z, '-q')
% z should be in the interval [-0.5 +0.5]
%
% Wb,May16,06

  N=0; if nargout>2, I=struct; end
  getopt ('INIT', varargin);
     if getopt({'nolog','-q'}), qflag=1;
     elseif getopt('-Q'), qflag=2; else qflag=0; end
     xi     = getopt('xi'   ); % return xi data relevant for NRGWilson.cc
     N      = getopt('N',N  ); % return xi data relevant for NRGWilson.cc
     ALambda= getopt('-AL'  ); % account for wavefunction renormalization A_Lambda
     Nmax   = getopt('Nmax',[]);
  getopt('check_error');       % in line with Oliveria (2005), Rok Zitko (2009)

  if numel(N)~=1 || isinf(N) || isnan(N)
     error('Wb:ERR','\n   ERR invalid N (%g)',N); end
  if abs(N-round(N))>0
     q=N; N=ceil(N);
     wblog('WRN','got N=%g -> N=%g',q,N);
  end

  if isempty(Nmax)
     betaMin=1E-16;
     Nmax = ceil( 2* abs(log(betaMin) / log(Lambda)) );
  else
     betaMin=0;
     if Nmax<0, Nmax = 2*N; end
  end

  if Gamma<0, error('Wb:ERR','invalid Gamma (%g)',Gamma); end

  if z==0 && ~ALambda
     ii=0:N-3; Delta=1/Lambda;
     xx=(1-Delta.^(ii+1))./sqrt((1-Delta.^(2*ii+1)).*(1-Delta.^(2*ii+3)));

     if xi
          b = [ sqrt(2*Gamma/pi), (1+Delta)/2                  * xx ];
     else b = [ sqrt(2*Gamma/pi), (1+Delta)/2*(Delta.^(ii/2)) .* xx ]; end
     a = zeros(1,length(b)+1);

     return
  elseif N>Nmax
     wblog('WRN','N too large; set N=%g -> %g (z=%g)', N, Nmax, z);
     N=Nmax;
  end

% NB! ensure continuous transition at z->0: z => 1-z
% otherwise for z<<1, small interval appears close to band edge!
% Wb,Feb17,11
  z=mod(1-z,1); % = mod(-z,1) :: z \el [0,1[ % tags: ZFLIP

  iz=0:Nmax; if z, iz=[0,iz+z]; end

  Edisc = (1/Lambda).^iz;
  n=length(Edisc)-1; % number of intervalls

% initialize Hamiltonian
  if ALambda % see pdf- and hand-notes (Zitko et al, 2009)
   % NB! there appears to be systematic error / offset for Lambda<2
   % e.g. standard Friedel sum rule for SIAM well-obeyed for Lambda>=2,
   % but A(0)=~0.98 becomes < 1 for Lambda=1.7 (!!) the latter case
   % seems to work better with ALambda turned off! (systematic error
   % compensates overshoot?) -- Wb,Jun16,09

   % NB! z=0 has narrow level at band-edge
   % -> actually prefer z->1 by default (smoother ff couplings).

     d0=((1-1/Lambda)/log(Lambda)) * Edisc(1:end-1);
     if z==0
        % use z=1 in expression below! Wb,Apr20,11 (hint Irek Weymann)
          d0(1)=      (1-1/Lambda  )/log(Lambda);
     else d0(1)=(1-z)+(1-1/Lambda^z)/log(Lambda);
     end

   % if z, d0(1)=(1-1/Lambda^z)/(z*log(Lambda));
   % else  d0(1)=1; end % Zitko (2009), see handnotes 
  else
     d0=0.5*(Edisc(1:end-1)+ Edisc(2:end)); % diagonal entries
  end

% coupling to impurity (Gamma included below) - tags: GAMMA
  d1 = sqrt(abs(diff(Edisc))); % * sqrt(Gamma/pi);

  hdd=sparse(diag(d0)); Z=sparse(n,n);

  H0=[ 0     d1    d1
       d1'  +hdd   Z
       d1'   Z    -hdd ];

  D=size(H0,1);    % 2*ndim + 1
  niter=D;         % max number of iterations

  if ALambda
     s='ALambda; ';
     global param; param.ALambda='using shifted discrete energies';
   % param.Eb=d0;
  else s=''; end

  if qflag % if nargout<2
     if qflag==1, wblog('<i>',...             % tags: ZFLIP
     'Oliveira (%sz=%g; Lambda=%g; Gamma=%g)',s,1-z,Lambda,Gamma); end
  else
     wblog('<i>','Oliveira (%sz=%g)',s,1-z); % tags: ZFLIP
     wblog('Lambda  = %g', Lambda);
     wblog('Gamma   = %g', Gamma );
     wblog('ndim    = %d -> %d', n, D);
  end

% starting lanczos tridiagonalization
  U=zeros(D,1); U(1)=1; % start vector

  alpha=zeros(1,niter); beta=zeros(1,niter);

  for i=1:niter
     v=H0*U(:,i);
     alpha(i)=U(:,i)'*v;

     v=v-U*(U'*v); v=v-U*(U'*v); % twice for numerical reasons
     beta(i)=norm(v);

     if beta(i)<=betaMin
        if ~qflag, wblog('<i>',...
        'Lanczos converged at iteration %d/%d (%.3g).',i,D,beta(i)); end
        break
     elseif i>=N, break;
     else U(:,i+1)=v/beta(i); end
  end

  alpha=alpha(1:i); beta=beta(1:i-1);

% NB! U(:,1) corresponds to normalized first column in H0 (!)
% apply missing factor to beta(1) - tags: GAMMA
  beta(1)=beta(1)*sqrt(Gamma/pi); % this also allows Gamma==0

  if xi
       a=alpha; b=beta .* [1, Lambda.^(((1:length(beta)-1)-1)/2)];
  else a=alpha; b=beta; end

% k=max(find(abs(beta)>1E-14));
  k=length(beta);
% NB! ignore ff(1) since it carries Gamma!
  br=abs(beta(3:k)./beta(2:k-1)); e=max(br);
  if e>=1
     % kx=max(find(br>0.9))+2; % shall descrease at least somewhat!
       kx=find(br==max(br),1);
  elseif beta(2)>beta(1), kx=2; else kx=1; end

  if kx>3
     s=sprintf(['got increasing couplings: ' ... 
       'ff(%g)=%.3g @ %.3g (z=%.3g)'],kx+1,beta(kx+1),br(kx),z);
     if e>2, error('Wb:ERR',['\n   ERR ' s]);
     else wblog(iff(kx<5,'WRN','ERR'),s); end
  end

  if nargout>2,
     I=add2struct('-',...
     Lambda,z,Gamma,Edisc,d0,alpha,beta,H0,U,kx,ALambda);

     if ALambda
      % see tst_oliveira_EScale.m // Wb,Apr26,11
        I.ALfac = Lambda^z*(Lambda-1)/log(Lambda);
     end
  end

  if qflag
     a=a(1:N); b=b(1:N-1);
     return
  end
  
% -------------------------------------------------------------------- %
% -------------------------------------------------------------------- %
% plot

ah=smaxis(2,2,'tag','LANCZOS');

  header(sprintf(...
  'NRG coefficients using Oliveira''s trick (z=%g, \\Gamma=%g, \\Lambda=%g)',...
   1-z, Gamma, Lambda)) % tags: ZFLIP

setax(ah(1,1));

  plot(abs(alpha), 'ko-'); hold on
  plot(abs(beta),  'ro-'); sms(4)

  plot(Edisc, 'b-'); % \DeltaE

  title(sprintf('\\Lambda=%g, \\Gamma=%g', Lambda, Gamma));
  legend('\alpha','\beta','\DeltaE'); legend boxoff
  set(gca,'YScale', 'log'); % ylim([1E-20 10]);
  axis tight

setax(ah(1,2)); % difference to Wilson chain

  beta0=zeros(size(beta));
  beta0(1)=sqrt(2*Gamma/pi);

  ii=2:length(beta);
  beta0(2:end) = ...
     (1+1/Lambda) ./ (2*Lambda.^((ii-2)/2)) .* ...
     (1-Lambda.^(-(ii-2)-1)) ./ (sqrt(1-Lambda.^(-2*(ii-2)-1)) .* ...
     sqrt(1-Lambda.^(-2*(ii-2)-3)));

  semilogy(abs(beta-beta0),'r'); hold on
  semilogy(abs(alpha),'k'); axis tight

  title('difference to Wilson''s \xi_{ i}')

setax(ah(2,1));

  T=diag(alpha)+diag(beta,1)+diag(beta,-1);
  
  mp0(H0-U*T*U','cmap','addrc'); 
  postext(0.1,0.95,'H-U*T*U'''); axis equal tight

setax(ah(2,2));

  mp0(U'*U-eye(size(U,2)),'cmap','addrc','diag',0);
  postext(0.1,0.95,'U*U''-eye()'); axis equal tight

  a=a(1:N); b=b(1:N-1);
return

% -------------------------------------------------------------------- %
% -------------------------------------------------------------------- %

