function se = SEntropy(r,varargin)
% Function SEntropy - calculate shanon entropy for given input vector
% Usage: se = SEntropy(r [,OPTS])
%
%    r are considered amplitudes, yet they should be positive.
%    before calculating SE, r is normalized.
%    SE is calculated using p_i = (r_i)^2
%
% Options
%
%    '-q'     quiet
%    'norm'   explicitely normalize to 1. (assumed otherwise)
%    'rho'    expect eigennvalues of reduced density matrix (SVD otherwise)
%
% Wb,Oct19,05  Wb,Aug05,06

  getopt('init',varargin);
     qflag  = getopt('-q');
     nflag  = getopt({'-n','norm'});
     gotsvd =~getopt({'-rho','rho'});
  getopt('check_error');

  if isempty(r), se=0; return; end
  if ~isvector(r)
     e=(r-r').^2; e=sum(e(:));
     if e/max(abs(r(:)))>1E-14
        wblog('ERR','invalid usage (size %s)',vec2str(size(r),'sep','x'));
        se=nan; if ~isbatch, keyboard, end, return
     end
     r=0.5*eig(r+r'); gotsvd=0;
     if ~isreal(r), wblog('WRN','eigenvalues not real !??'); end
  end
  if ~qflag && min(r)<-1E-14
     wblog('WRN','negative entries for shannon entropy ?? (%g)', min(r));
  end

  r=r(:); if nflag, r=r/norm(r); end

  if gotsvd, r=r.*conj(r); end
  r(find(r<=0))=[]; % avoid warnings by 0*log(0)

  if ~qflag
     e=sum(r); if abs(e-1)>1E-10
     wblog('WRN','rho_i not normalized !?? (1%+.3g)',e-1); end
  end

  se = -r'*log(r);

end

