function [TK,LK] = TKondo(varargin)
% function: TKondo
% usage   : [TK,LK] = TKondo
%
%    takes param or paras as input and calculates Kondo temperature
%    
%    in addition, given a certain logarithmic discretization specified
%    by param.Lambda, LK gives the chain length L needed to resolve TK
%
% Wb,Mar01,05

  global g_TKondo

  if ~isempty(g_TKondo)
     if ~isnan(g_TKondo)
     wblog(' * ','using global TK=%.3g', g_TKondo); end
     TK=g_TKondo; check_value(TK);
     return;
  end

% if(checkParam(0)) return; end

  getopt('init',varargin);
     lflag=getopt('log');
     wTK  =getopt('wtk',1);
  varargin=getopt('get_remaining'); narg=length(varargin);

  if narg, param=varargin{1};
     if narg>1 || ~isstruct(param), error('Wb:ERR','invalid usage'); end
  else global param
  end

  structexp(param);

  if exist('fJ','var') && ~exist('J','var'), J=fJ; end
  if exist('J','var') && ~exist('epsd','var')
     if J==0, TK=0; check_value(TK); return; end

   % TK=sqrt(2*J)*exp(-1/(2*J)); % Hewson, p.63, (3.47))

   % effective matrix element for single channel spin 1/2 due
   % to coupling: J/4 (Costi, PRL,2000)
   % TK=sqrt(J/2)*exp(-1/(J/2));

   % in Rosch et al. (PRB 2003) the matrix elements are J/2
   % same as in M.Sindel (thesis p. 36)
   % equation TK is offset from energy flow convergence regime
   % => take J->J/2 for definition of Kondo temperature
   % NB!  Rosch et al. (PRB 2003) footnote [33]
   %   TK =~ sqrt(J)*exp(-1/J)   - exact within 2-loop perturbation
   %   TK =~         exp(-1/J)   - exact within 1-loop perturbation

     if J>0
     TK=sqrt(J)*exp(-1/J); % using J = fJ = J_calc = 2 \nu JK = 2nJ
                           % ==> TK \equiv sqrt(2nJ) exp(-1/2nJ) (!!)
   % equivalently: using (2D nuJ) Sd.S0 = (nuJ) Sd.tau0 as Kondo Hamiltonian
   % => J_NRG = nuJ (see setupKondo*.m)
   % => still: TK =~ sqrt(nuJ)*exp(-1/nuJ) // Wb,Sep08,13
     else TK=nan; end

     check_value(TK);
     return
  end

% actually, max(TK) is obtained at symmetric point epsd=-U/2!
% if (epsd+U)<0, wblog('WRN','epsd+U<0 !?? - take |epsd+U|.'); end

  lastwarn(''); LK=[];

  if ~exist('epsd','var') || ~exist('Gamma','var') || ~exist('U','var')
     TK=nan; check_value(TK);
     return;
  end

  try
    e=epsd;  if length(e)>1, e=mean(e(:)); end
    g=Gamma; if length(g)>1, g=mean(abs(g(:))); end
    u=U;     if length(u)>1, u=mean(u(:)); end

    if g==0, TK=0; check_value(TK); return; end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % * Hewson, p.66
  %   jk =abs(e*(e+u))/(g*u); TK = sqrt(1/jk)*exp(-pi*jk/2);
  %   TKondo = sqrt(U*Gamma) * exp(pi*epsd*(epsd+U)/(U*Gamma));
  % * Thesis, M. Sindel, p.29
  %   NB! requires Gamma -> Gamma/2, otherwise TK too small!)
  % * Haldane, PRL, 6, 416 (1978) - Delta \equiv myGamma
  %   TK = sqrt(u*g)/2 * exp(-abs(pi*e*(e+u)/(u*g)));
  % * for symmetric Anderson model (see Merker12 + Costi) // Wb,Jun20,12
  %   TK = sqrt(u*g/2) exp(-pi(u/8g - g/2u)) \equiv 1/4chi0
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % TK = sqrt(u*g/2) * exp(pi*e*(e+u)/(2*u*g));
    if u==0, TK=0; else
       e=pi*e*(e+u)/(2*u*g); TK=nan; wTK=1; % default
       switch wTK

       % -------------------------------------------------------- %
       % NB! the prefactor of sqrt(u*g/2) is somewhat problematic
       % as it goes to zero for U->0 (resonant-level model) (!?!)
       % Wb,Mar28,13
       % -------------------------------------------------------- %
         case 1, TK=min(0.575,sqrt(u*g/2))*exp(e);
         % actually using 2*Gamma instead of Gamma (!?)
         % consistent with definition of Costi '94 (!!)
         % prefactor chosen to agree with TK_chi for large U // Wb,May08,13

         case 2, TK=sqrt(u*g/4)*exp(2*e);
         % Gamma->Gamma/2 (!) eg. see review by Glazman (?)
         % then TK is somewhat small (indicates already
         % converged Kondo regime?)

         case 3, TK=sqrt(u*g/(2*3.00))*exp(1.050*e);
         % find best overlap with magn. suszeptibility

         otherwise error('Wb:ERR','invalid usage');
       end

       if isinf(TK), TK=nan; end
    end
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    if isnan(TK) && u && e<0
       warning('Wb:ERR:TK','U*Gamma=%.3g -> TK is NaN', g*u);
    elseif TK>10
     % warning('Wb:ERR:TK','epsd=%.3g U -> TK=%.4g !??', e/u,TK);
       TK=nan;
    elseif TK<0
       warning('Wb:ERR:TK','Gamma=%.3g U -> TK=%g', g/u, TK);
       TK=nan;
    end
  catch
    l=lasterror; dispstack(l);
    warning('Wb:ERR:TK', strrep(l.message, char(10), '\n'));
  end

  if ~exist('TK','var'), TK=nan; check_value(TK); return; end

  if isnan(TK), check_value(TK); return; end
  if ~isempty(lastwarn)
     wblog('WRN','Failed to evaluate expression for TK');
     wblog('CST',sprintf('MSG %s', lastwarn)); % keyboard
     return
  end

  if nargout<2; check_value(TK); return, end

  if exist('Lambda')
  LK = -2 * log(TK)/log(Lambda); end

  if lflag, wblogx(1,...
  '<i> TKondo = %e (LK=%.1f)', TK,LK); end

end

% -------------------------------------------------------------------- %
% -------------------------------------------------------------------- %
 
function check_value(TK)

   if isnan(TK) || TK>10 || TK<=0
      wblog('WRN','TK=%g',TK);
   end

end

% -------------------------------------------------------------------- %
% -------------------------------------------------------------------- %
 
