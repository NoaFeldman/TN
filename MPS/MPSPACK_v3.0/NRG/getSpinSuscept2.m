function [chi,Iout]=getSpinSuscept2(Simp,varargin)
% function [chi,I]=getSpinSuscept2(Simp [,nrg,opts])
%
%    calculate mixed spin susceptibility derived from Matsubara
%    expression for <Simp||Stot> = <Simp*Stot>. Here Simp is the spin
%    operator at the impurity represented in the local basis of A0.
%
% Wb,Jun19,13

% adapted from getSpinSuscept2_iter.m
% see also rnrg_tst_chi2.m

  getopt('init',varargin);
     vflag=getopt('-v'); if ~vflag
     qflag=getopt('-q'); else qflag=0; end
     TT=getopt('T',[]);
     iS=getopt('iS',{}); if ~isempty(iS), iS={iS}; end
  varargin=getopt('get_remaining');

  if numel(varargin) && ischar(varargin{1})
       nrg=varargin{1}; varargin=varargin(2:end);
  else nrg=[getenv('LMA') '/NRG/NRG']; end

  if numel(varargin) || ~ischar(nrg)
     helpthis, if nargin || nargout
     error('Wb:ERR','invalid usage'), end, return
  end

  f=sprintf('%s_info.mat',nrg);
  if ~exist(f,'file')
     error('Wb:ERR','\n   ERR invalid NRG data (%s)',nrg);
  end
  Inrg=load(f);

  if isempty(TT)
     mt=4;
     TT=(250*Inrg.EScale(end))*Inrg.Lambda.^((0:mt-1)/mt);
  end

  nT=numel(TT); rw=zeros(Inrg.N,nT);

  i=find(TT<=0); n=numel(i);
  if n
     if n>1, error('Wb:ERR',...
       '\n   ERR invalid T (nan may only be set once)'); 
     end
     rw(:,i)=repmat(Inrg.rhoNorm(:,2),1,n);
     TT(i)=Inrg.rhoT;
  end

  if n<nT
     if vflag, o={'-v'}; else o={}; end
     if ~qflag, wblog(' * ',...
       'get FDM weights for %g temperatures',numel(TT));
     end

     [rw,TT]=getWeightsNRG(nrg,'T',TT,o{:});
     if n
        e=norm(rw(:,i)-repmat(Inrg.rhoNorm(:,2),1,n));
        if norm(TT(i)-Inrg.rhoT)~=0 || e>1E-12
           error('Wb:ERR','\n   ERR FDM weigth inconsistency (@ %.3g)',e); 
        end
     end
  end

  i=find(rw(end,:)>0.01);
  if numel(i)
     j=find(rw(end,:)==max(rw(end,:)));
     if numel(j)>1, j=j(find(TT(j)==min(TT(j)),1)); end
     wblog('WRN','got too small a temperature (T=%.3g @ %.3g) !??',...
     TT(j),rw(end,j));
  end

  N=size(rw,1); chi=zeros(nT,N);

  if ~qflag, wblog(' * ',...
    'calculate spin susceptibility <Simp||Stot> ...',numel(TT));
  end

% NB! use k as outer loop, so SKK is calculated once
% for all temperatures!
  for k=1:N
     if vflag
        fprintf(1,'\r   k=%2g/%g (T=%.3g .. %.3g; %g values) ... \r',...
        k,N,min(TT),max(TT),nT);
     end
     q=load(sprintf('%s_%02g.mat',nrg,k-1),'AK','AT','HT','ES');

     if k>1
        if ~isempty(q.AT.data)
           SDD=contractQS(q.AT,'*12',{q.AT,SKK},'-L');
        end
        if k<N
           SKK=contractQS(q.AK,'*12',{q.AK,SKK},'-L');
        end
     else
        if ~isempty(q.AT.data), error('Wb:ERR',['\n   ERR got ' ... 
           'NRG truncation already with first Wilson site !??']);
        end
        d=dim(QSpace(q.AK),'-a');
        if d(end,1)>1
             SKK=contractQS(q.AK,'*12',{q.AK,Simp,'-op:L00'}); kimp='L';
        else SKK=contractQS(q.AK,'*23',{q.AK,Simp,'-op:s00'}); kimp='s';
        end

        if isempty(SKK.data)
           error('Wb:ERR','\n   ERR invalid Simp contraction (empty)');
        end
     end

     if isempty(q.AT.data), continue; end

     for it=1:nT
        if vflag
         % fprintf(1,'\r   k=%2g/%g: it=%g/%g ...  \r',k,N,it,nT);
        end

        if ~isempty(SDD.data)
         % in case that the discarded space is very small
         % (e.g. at first iteration with truncation, it can 
         % happen tath SDD is empty!) // Wb,Jun29,13
           chi(it,k) = (rw(k,it)/q.ES) * ...
           getSpinSuscept2_iter(SDD,q.HT, TT(it)/q.ES, iS{:});
        end
     end
     if vflag, fprintf(1,'\r%50s\r',''); end
  end

  if nargout>1
     Iout=add2struct('-',chi,rw,TT,kimp);
  end

  chi=reshape(sum(chi,2),size(TT));

end

% -------------------------------------------------------------------- %
% -------------------------------------------------------------------- %

function chi=getSpinSuscept2_iter(Simp,H,T,iS)
% function chi=getSpinSuscept2_iter(Simp,H,T,iS)
%
%    calculate (contribution to) mixed spin susceptibility
%    <Simp||Stot> = <Simp*Stot>. Here Simp is the matrix elements
%    of S_imp represented in the effective basis of given NRG
%    iteration, H contains the energy of each state, which
%    is required for the construction of the thermal density
%    matrix at temperature T.
%
%    For the construction of Stot, the position
%    of the SU(2) spin symmetry in info.qtype is required.
%    if iS is not specified, this routine tries to
%    determine iS from the structure of Simp.
%
% Wb,Jun19,13

% adapted from rnrg_tst_chi2.m (for further tests, see there).

  if nargin==1 && isa(Simp,'QSpace')
   % for testing purposes: explicitly construct and return Stot
     chi=get_Stot(Simp); return
  elseif nargin<3 || nargin>4
     helpthis, if nargin || nargout
     error('Wb:ERR','invalid usage'), end, return
  end

% Simp_=Simp;
  Simp=struct(Simp);

  if nargin<4
     if numel(Simp.Q)<3, error('Wb:ERR',[...
       '\n   ERR cannot determine iS from abelian Simp (Sz only)' ...
       '\n   ERR => must specify iS upon input']);
     end

   % look for symmetry labels with largest variations
     q=std((Simp.Q{1}-Simp.Q{2}).^2,[],1); i=find(q); n=numel(i);
     if n>1
        s=sprintf('got %g symmetries with varying q-labels (%g)',n);
        if numel(find(q>0.99*max(q)))>1
             error('Wb:ERR',['\n   ',s]);
        else wblog('WRN',s); end
     end
     if max(q)>1E-12
        iS=find(q==max(q));
     else
        iS=(find(Simp.Q{3}(1,:)==2));
     end
     if numel(iS)~=1, iS
        error('Wb:ERR','\n   ERR failed to determine iS');
     end
  else
     if numel(iS)~=1, iS
        error('Wb:ERR','\n   ERR invalid iS');
     end
  end

  if numel(Simp.Q)>2
     if norm(diff(Simp.Q{3},[],1))>1E-8 || norm(Simp.Q{3}(:,iS)-2)>1E-8
        error('Wb:ERR',['\n   ' ... 
       'ERR invalid Simp (expecting iROP here; %g)!'],iS);
     end
  end

% -------------------------------------------------------------------- %
% build Stot-operator *implicitly* from Simp
%
% (1) Simp already has all the correct CGC spaces!
% (2) note that Stot is *diagonal* in given effective state space!
%     that is, Q{1}==Q{2} must be equal throughout! (this is
%     justified by the fact S_+ only affects the Sz label of
%     *given* multiplet!)
%
%  => can reduce Simp to block-diagonal
%     (only these enter in the final contraction).
% -------------------------------------------------------------------- %

  i=find(sum(abs(Simp.Q{1}-Simp.Q{2}),2)<1E-8);
  Simp=getsub(QSpace(Simp),i);

  dS=dim(Simp,'-a'); dS(:,end+1:3)=1;
  dop=dS(end);

  aflag=iff(diff(dS(:,3)),0,1); % whether got abelian Spin

  if ~isa(H,'QSpace'), H=QSpace(H); end
  if isdiag(H,'-d')~=2, error('Wb:ERR',...
    '\n   ERR invalid Hamiltonian (expected in diagonal format!)');
  end

  [i1,i2,I]=matchIndex(Simp.Q{1},H.Q{1});
  if ~isempty(I.ix1), error('Wb:ERR',...
     '\n   ERR Simp appears incomplete/incompatible w.r.t. given H!');
  end

  if T<=0, error('Wb:ERR','\n   ERR invalid temperature (%g)',T); end
  beta=1/T;

  R=getrhoQS(H,beta);
  X=contractQS(R,2,Simp,1);

% contract Stot by modifying X=(R*Simp):
% since Stot is (block-)diagonal, only plain factors arise
% which can be contracted onto X.data{:}, while already
% also tracing over the latter space.

  nq=size(X.info.cgs,2);

  for i=1:numel(X.data), S=X.Q{1}(i,iS)/2;
     if aflag, a=S;
     else
      % NB! SS for local space has minus sign;
      % hence, for consistency, also use reversed sign here!
        a=-sqrt(S*(S+1));
     end
     m=1; % multiplicity (i.e. z-space dimension)

   % contract CGC to obtain correct weight for Stot!
     for j=1:nq, q=X.info.cgs(i,j);
        if isfield(q,'S')
           if isempty(q.S)
              if ~isequal(q.data,1), q, error('Wb:ERR',...
                '\n   ERR unexpected abelian CGC space'); end
              continue
           end
        else % Wb,Sep21,13
           if ~isequal(q,{speye(size(q{1},1))}), q, error('Wb:ERR',...
             '\n   ERR unexpected abelian CGC space'); end
           continue
        end

        m=m*q.S(2); q=reshape(full(q.data),q.S);
        q=contract(q,q,[1 3],[1 3]);
        d=diag(q); if std(d)>1E-12
           wblog('ERR','unexpected CGC space (%.3g)',std(d)); end
      % NB! for standard CGC's, d=1 (!)
      % already include factor from the trace of X with Stot
      % a=a/sqrt(mean(d)) * mean(d); // see rnrg_tst_chi2.m
        a=a*sqrt(mean(d));
     end
     if abs(a)<1E-8 && S~=0 || m<1
        error('Wb:ERR','\n   ERR invalid CGC fac (%.3g; %g)',a,m);
     end

     X.data{i}=(m*a)*trace(X.data{i});
  end

  chi=(beta/dop)*sum(cat(1,X.data{:})); % = beta * trace( R * Simp * Stot )
  % dop=dim(op) accounts for IROP factor

end

% ==================================================================== %
% for archiving purposes (the following function is not called)
% -------------------------------------------------------------------- %
% build Stot-operator from Simp at last iteration
% NB! already has the correct CGC spaces!

% adapted from rnrg_tst_chi2.m
% where norm(Stot - Stot_exact)=0, indeed; there Stot_exact
% was explicitly contstructed as the iteratively calculated 
% matrix elements of Stot = sum_i Stot_i.

function Stot=get_Stot(Simp)

  Stot=Simp(end,1);
  i=find(sum(abs(Stot.Q{1}-Stot.Q{2}),2)<1E-8);
  Stot=getsub(Stot,i);

  iS=2; Stot=struct(Stot);
  for i=1:numel(Stot.data), S=Stot.Q{1}(i,iS)/2;
     q=Stot.data{i}; s=min(size(q));
   % NB! SS for local space has minus sign;
   % hence, for consistency, also use reversed sign here!
     x=-sqrt(S*(S+1));

   % contract CGC to obtain correct weight
     for j=1:size(Stot.info.cgs,2)
        q=Stot.info.cgs(i,j); q=reshape(full(q.data),q.S);
        q=contract(q,q,[1 3],[1 3]); d=diag(q);
        if std(d)>1E-12, wblog('ERR','unexpected CGC space'); end
        x=x/sqrt(mean(d));
     end
     Stot.data{i}=diag(repmat(x,1,s));
  end
  Stot=skipzeros(QSpace(Stot));

end

% ==================================================================== %
% -------------------------------------------------------------------- %

