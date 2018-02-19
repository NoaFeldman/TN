function [se,rr,I]=SEntropy(R,varargin)
% function [se,rr,I]=SEntropy(Rho [,eps,opts])
%
%   calculate shannon entropy where the reduced density matrix R
%   also may contain non-abelian symmetries. In the latter case,
%   the output weights rr are already also weighted by their 
%   corresponding multiplet degeneracies.
%
% Options
%
%    eps      check trace(Rho)=1 within numerical accuracy of eps (1E-12)
%   '-xm'     expand multiplet degeneracy in returned rr
%   'beta',.. assume R is a Hamiltonian => apply beta=1/T and exponentiate
%
% Wb,Jan24,11

  if nargin<1 || numel(R)~=1
     eval(['help ' mfilename]);
     if nargin || nargout, error('Wb:ERR','invalid usage'), end, return
  end

  getopt('init',varargin);
     xmflag=getopt('-xm');  % Wb,May09,14
     beta=getopt('beta',0); % Wb,May27,16
  eps=getopt('get_last',1E-12);

  if numel(R.Q)~=2 || ~isequal(R.Q{1},R.Q{2})
  error('Wb:ERR','invalid density matrix'); end

  if beta>0, nd=numel(R.data);
     if isdiag(R)>1 && nargout<2   % Wb,May27,16
      % no need to diagonalize using eigQS() below!
        drs=get_qdim(R,2);
        for i=1:nd
           R.data{i}=exp(-beta*R.data{i});
           drs(i,2)=sum(R.data{i});
        end
        R=R/sum(prod(drs,2));
        for i=1:nd, r=R.data{i}; j=find(r);
           if ~isempty(j), rj=r(j);
              drs(i,3)=-rj*log(rj)';
           end
        end
        se=drs(:,1)'*drs(:,3); % keyboard
        return
     end
     for i=1:nd, r=R.data{i}; s=size(r);
        if any(s==1),      r=exp (-beta*r);
        elseif diff(s)==0, r=expm(-beta*r);
        else error('Wb:ERR','\n   ERR invalid operator H'); 
        end
        R.data{i}=r;
     end
   % see normalization (incl. CGC factors!) below
  end

  if isdiag(R)>1, R=diag(R); end
  if beta>0, R=R/trace(R); end

  check_trace(trace(R),eps,nargin>1);

  [rr,I]=eigQS(R);

     if size(rr,2)==2, rr=rr(:,1).*rr(:,2); end % Wb,Aug25,15
     i=find(rr<=0); e=norm(rr(i));
     if e>1E-12, error('Wb:ERR','invalid density matrix (%g)',e); end

% if isempty(R.info)
  if isAbelian(R)
   % plain setting for abelian symmetries, i.e. not CGS space present
     rr=sort(rr,'descend');
     se=SEntropy(rr,'-rho'); return
  end

% account for dimension in CGC space
% if any(diff(rr)<0), error('Wb:ERR','invalid usage'); end
  if numel(rr)~=sum(I.DB(:,1)), error('Wb:ERR','dimension mismatch'); end

  if size(I.DB,2)>1
     % total dimension incl. CGC / total multiplet dimension = CGC dimension
       dz=I.DB(:,2)./I.DB(:,1); % combined dimension of CGC space
  else dz=ones(size(I.DB)); end

% rr=mat2cell(rr,1,dz); NB! rr above is sorted!
  rr=I.EK.data;

  I.rr=rr;
  I.dz=dz; se=0; n=numel(dz);

  for k=1:n
     r=rr{k}; r(find(r<=1E-16))=[];
     se = se - dz(k)*(r*log(r)');
     if ~xmflag
        rr{k}=dz(k)*rr{k}; % default
     else
      % expand eigenvalues to level of state space which makes it
      % easier to compare data with different symmetry settings!
      % eg. see NRG/getDiscWeightNRG.m // Wb,May11,14
        rr{k}=repmat(rr{k},1,dz(k));
     end
  end

  [rr,is]=sort(cat(2,rr{:}),'descend');

% keyboard
  check_trace(sum(rr),eps,nargin>1);

end

function check_trace(t,eps,fflag)

  if norm(1-t)>eps
     s=sprintf('invalid density matrix (tr(rho)=%.3g @ %.3g)',t,1-t);
     if fflag || norm(t-1)>0.5, error('Wb:ERR',s);
     else 
        if norm(t-1)<1E-4, wblog(1,'WRN',s); else wblog(1,'ERR',s); end
     end
  end

end

