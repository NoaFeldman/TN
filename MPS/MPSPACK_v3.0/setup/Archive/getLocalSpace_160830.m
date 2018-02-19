function varargout=getLocalSpace(model,varargin)
% function [FF,..,Io]=getLocalSpace(model [,'sym1,sym2,...',varargin])
%
%   build local model state space as specified by means of the typical
%   associated operators such as spin (S), annihilation (F), or
%   charge parity operator (Z).
%
%   The residual info structure Io contains further operators
%   (if applicable), such as the identity operator (E), the spinor for
%   particle-hole symmetry (C3), its equivalent to the Casimir operator
%   S^2 (C2), or Q2 := (N-1)^2, with N the total particle number summed
%   over all channels, as required for isotropic Coulomb interaction.
%
% Models available
%
%    [F,Z,S,I]=getLocalSpace('FermionS',sym [,opts]);   % spinful fermions
%    [F,Z,  I]=getLocalSpace('Fermion', sym [,opts]);   % spinless fermions
%    [S,    I]=getLocalSpace('Spin', S      [,opts]);   % spin-S operators
%    [S,    I]=getLocalSpace('SUN',  N      [,opts]);   % SU(N) site [so far: N=3]
%
% Options
%
%   'NC',..  number of channels (fermionic systems only)
%   '-A',..  use abelian symmetry (in 'Spin' mode only)
%   '-v'     verbose mode
%
% Symmetries for sym (single string, separated by commas)
%
%   'Acharge'      abelian total charge;        Acharge(:)  *)
%   'SU2charge'    SU(2) total particle-hole;   SU2charge(:) *)
%   'Aspin'        abelian total spin
%   'SU2spin'      SU(2) total spin
%   'SUNchannel'   SU(N) channel symmetry
%
% *) sym(:) indicates to use given symmetry <sym> for each
%    of the NC channels individually.
%
% Examples
%
%  % three channels with particle-hole SU(2) in each channel and total spin SU(2)
%    [FF,Z,SS,IS]=getLocalSpace('FermionS','SU2charge(:),SU2spin','NC',3,'-v');
%  % three spinless channels with SU(3) channel symmetry
%    [FF,Z,IS]=getLocalSpace('Fermion','SUNchannel','NC',3,'-v');
%  % two channels with abelian charge, SU(2) spin and SU(2) channel
%    [FF,Z,SS,IS]=getLocalSpace('FermionS','Acharge,SU2spin,SUNchannel','NC',2,'-v','-Y');
%  % single spin-S site
%    [S,IS]=getLocalSpace('Spin',1,'-v'); 
%
% Wb,Jul09,11 ; Wb,Jul30,12

% ---------------------------------------------------------------------- %
% for testing purposes
% * two channels with particle-hole and spin SU(2) in each channel
%      [FF,Z,SS,IS]=getLocalSpace(...
%      'FermionS','SU2charge(:),SU2spin(:)','NC',2,'-v');
% ---------------------------------------------------------------------- %
 
% fprintf(1,'\n'); eval(['help ' mfilename]);
  if nargin<1 || ~ischar(model)
     helpthis, if nargin || nargout
     error('Wb:ERR','invalid usage'), end, return
  end

  if nargin>1 && isequal(varargin{end},'-fixA')
       fixA=1; varargin=varargin(1:end-1);
  else fixA=0; end

  if ~isempty(regexp(model,'^SU\d$'))
     varargin=[{str2num(model(3:end))}, varargin];
     model='SUN';
  end

  switch model
     case 'FermionS'
        n=4; varargout=cell(1,n); % {F,Z,S,Io}
        [varargout{:}]=getLocalSpace_SpinfullFermions(varargin{:});

     case 'Fermion'
        n=3; varargout=cell(1,n); % {F,Z,Io}
        [varargout{:}]=getLocalSpace_SpinlessFermions(varargin{:});

     case 'Spin'
        n=2; varargout=cell(1,n); % {S,Io}
        if nargin>1 && ischar(varargin{1}) && regexp(varargin{1},'^SU\d')
             [varargout{:}]=getLocalSpace_SpinSUN(varargin{:});
        else [varargout{:}]=getLocalSpace_Spin(varargin{:});
        end

     case 'SUN'
        n=2; varargout=cell(1,n); % {S,Io}
        [varargout{:}]=getLocalSpace_SUN(varargin{:});

     otherwise error('Wb:ERR','\n   ERR invalid model ''%s''',model);
  end

      if nargout<2, varargout=varargout(1);
  elseif nargout<n, varargout=varargout([1:nargout-1,end]);
  elseif nargout>n, error('Wb:ERR',['\n  ERR invalid usage: ' ... 
     'asking for too many output arguments (%g/%g)'],nargout,n); end

% strip all info/cgr data from all-abelian QSpace operators
% Wb,Apr09,13
  if fixA
  for io=1:numel(varargout), q=varargout{io}; m=0;
     if isa(q,'QSpace')
        if any(isAbelian(q)==2)
           q=fixAbelian(q); m=m+1;
        end
     elseif isstruct(q), ff=fieldnames(q);
        for j=1:numel(ff), s=getfield(q,ff{j});
           if isa(s,'QSpace') && any(isAbelian(s)==2)
              q=setfield(q,ff{j},fixAbelian(s)); m=m+1;
           end
        end
     end
     if m, varargout{io}=q; end
  end
  end

end

% -------------------------------------------------------------------- %
% -------------------------------------------------------------------- %
% centralize initialization of qfac and jmap!
% Wb,Apr04,14

function q=init_qstruct(istr,sym,N)

   if nargin<2 || nargin>3 || ~ischar(sym) || ~ischar(istr)
      error('Wb:ERR','invalid usage (init_qstruct)'),
   end

   q.info=istr;
   q.type=[]; % see below

   q.Sp=SymOp; % to be set in caller
   q.Sz=SymOp;

   q.qfac=[]; % factor(s) applied on z-labels
   q.jmap=[]; % subsequent map for J-labels, ie. Jout => J_out*qfac

   switch sym
    case 'SU'
       if nargin<3, error('Wb:ERR',...
         '\n   ERR invalid usage (missing N)'); end
       if numel(N)~=1 || N<2 || N>10, error('Wb:ERR',...
         '\n   ERR invalid symmetry ''%s(%g)''',sym,r); end

       r=N-1; if r>1
        % SU(3) => r=2: q.jmap=[1 0  ; -1/2 1/2              ]';
        % SU(4) => r=3: q.jmap=[1 0 0; -1/2 1/2 0; 0 -1/3 1/3]';
          q.jmap=[ diag(1./(1:r)) - diag(1./(2:r),-1) ]';
       end

       q.type=sprintf('%s%g',sym,N);
       q.qfac=2; % also required for N=3

    case 'Sp' % tags: Sp4 Sp6 Sp8 // SP4 SP6 SP8

       if nargin<3, error('Wb:ERR',...
         '\n   ERR invalid usage (missing N)'); end
       if numel(N)~=1 || N<2 || N>20 || mod(N,2), error('Wb:ERR',...
         '\n   ERR invalid symmetry ''%s(2*%g)''',sym,N/2); end
       if N==2, wblog('WRN',...  % NB! Sp2 = SU(2)
         'got Sp(%g) which is equivalent to SU(2)',N); end

       r=N/2;

     % if r==3 % tags: Sp6 SP6
     %  % see MATH//tex-notes // SPN_LIE
     %  % see also clebsch.cc::findMaxWeight()!
     %    q.jmap=[6 0 0; -3 3 0; 0 -2 2]'/6; % NB! same as for SU3 above!
     % end

       if r>1
          q.jmap=[ diag(1./(1:r)) - diag(1./(2:r),-1) ]';
       end

       q.type=sprintf('%s%g',sym,N);

    otherwise
      
       if nargin>2, error('Wb:ERR',['\n   ERR invalid usage ' ... 
         '(got 3rd argument, having ''%s'')'],sym); end
       q.type=sym;
   end

end

% -------------------------------------------------------------------- %
% -------------------------------------------------------------------- %

function [F,Z,S,Iout]=getLocalSpace_SpinfullFermions(Sym_,varargin)
% function [F,Z,S,Io]=getLocalSpace_SpinfullFermions(Sym,varargin)
% supported symmetries
%
%   'Acharge'      abelian total charge
%   'SU2charge'    SU(2) total particle-hole
%   'Aspin'        abelian total spin
%   'SU2spin'      SU(2) total spin
%   'SUNchannel'   SU(N) channel symmetry
%   'SpNchannel'   Sp(N) channel symmetry
%
% Example: 'Acharge,SU2spin'
% tags: symplectic symmetry, Sp2n, Sp4, Sp6, SpNchannel
% Wb,Jul09,11

  if nargin<1 || ~ischar(Sym_)
     helpthis, if nargin || nargout
     error('Wb:ERR','invalid usage'), end, return
  end

  getopt('INIT',varargin);
     NC   =getopt('NC',1);
     Yops =getopt('-Y');

     if getopt('-v'), vflag=1; oc={'-v'};
     elseif getopt('-V'), vflag=2; oc={'-V'};
     else vflag=0; oc={}; end

  getopt('check_error');

  istr=sprintf('fermionic %g-level system',NC);

  f2=sparse([0 1; 0 0]); % => [1,0] = empty state
  z2=f2*f2'-f2'*f2; n2=f2'*f2; e2=speye(2);

% build state space: [ 1up,1down; 2up,2down; ...; NCup,NCdown ]
% F-operators
  FF=SymOp(NC,2); NN=SymOp(NC,2); SS=SymOp(NC,3);
  CP=SymOp(NC,1); CZ=SymOp(NC,1); n2=2*NC;

% NB! use -sym flag, such that also in abelian case S^2 = S.S (!)
  sx=spinmat(2,'-sp','-sym');
  sx=sx([1 3 2]); % [ Sp, Sm, Sz ] order

  for i=1:NC
     for j=1:2, k=2*(i-1)+j;
        q=[ repmat({z2},1,k-1), {f2}, repmat({e2},1,n2-k) ];
        FF(i,j)=SymOp(getOpName('fops',i,j),mkron(q{:}));
        NN(i,j)=FF(i,j)'*FF(i,j); NN(i,j).istr=sprintf('N%g',i);
     end

     k=2*i-1;
        q=[ repmat({e2},1,k-1), {mkron(f2',f2')}, repmat({e2},1,n2-k-1) ];
        CP(i)=SymOp(sprintf('C%g(+)',i),mkron(q{:})); q{k}=0.5*comm(q{k},q{k}');
        CZ(i)=SymOp(sprintf('C%g(z)',i),mkron(q{:}));

     for j=1:3
        q=[ repmat({e2},1,k-1), {blkdiag(0,sx{j},0)}, repmat({e2},1,n2-k-1) ];
        SS(i,j)=SymOp(getOpName('spin',i,j),mkron(q{:}));
     end
  end

  if Yops
     CC=SymOp(n2,n2); Id=SymOp('Id',speye(size(FF(1).op)));
     for i=1:n2
        for j=1:n2
           CC(i,j)=FF(i)'*FF(j);
           if i==j
            % NB! ensure traceless and orthogonal operator space!
            % (set PHYS notes: operator basis to Schrieffer-Wolff)
            % Wb,Sep28,12
              CC(i,j)=CC(i,j)-0.5*Id;
           end
        end
     end
  else
     CC=SymOp(0,0);
  end

  SOP=[]; D=dim(FF); phflag=0; chflag=0;

  Sym=strread(Sym_,'%s','whitespace',';, ')'; clear q

  for s=Sym, s=s{1};
  % change SU<N>channel to generall 'SUNchannel' while checking N==NC
    n=str2num(regexprep(s,'SU(\d+)channel','$1'));
    if ~isempty(n)
       if ~isequal(n,NC), error('Wb:ERR',...
         '\n   ERR invalid symmetry %s (having NC=%g)',Sym_,NC); 
       end
       s='SUNchannel';
    end

    n=str2num(regexprep(s,'Sp(\d+)channel','$1'));
    if ~isempty(n)
       if ~isequal(n,NC), error('Wb:ERR',...
         '\n   ERR invalid symmetry %s (having NC=%g)',Sym{i},NC); 
       end
       s='SpNchannel';
    end

    if isequal(s(end-2:end),'(:)') % apply symmetry on each individual channel
         cflag=1; if NC==1, s=s(1:end-3); end % irrelevant
    else cflag=0; end

    switch s
      case 'Acharge'    % abelian total charge
          q=init_qstruct('total charge U(1)','A');
        % q.Sp=SymOp;
          q.Sz=2*sum(CZ);
          if ~isempty(SOP), SOP(end+1)=q; else SOP=q; end

      case 'Acharge(:)'  % abelian charge for each channel
          for i=1:NC
             q=init_qstruct(sprintf('charge U(1) [%d/%d]',i,NC),'A');
           % q.Sp=SymOp;
             q.Sz=2*CZ(i);
             if ~isempty(SOP), SOP(end+1)=q; else SOP=q; end
          end

      case 'Aspin'      % abelian total spin

          q=init_qstruct('total spin U(1)','A');
        % q.Sp=SymOp;
          q.Sz=2*sum(SS(:,3));
          if ~isempty(SOP), SOP(end+1)=q; else SOP=q; end

      case 'SU2charge'  % SU(2) total particle-hole

          q=init_qstruct('total particle-hole SU(2)','SU',2);
          q.Sp=sum(CP);
          q.Sz=sum(CZ);
          if ~isempty(SOP), SOP(end+1)=q; else SOP=q; end
          phflag=1;

      case 'SU2charge(:)'  % SU(2) particle-hole for each channel

          for i=1:NC
             q=init_qstruct(sprintf('particle-hole SU(2) [%d/%d]',i,NC),'SU',2);
             q.Sp=CP(i);
             q.Sz=CZ(i);
             if ~isempty(SOP), SOP(end+1)=q; else SOP=q; end
          end
          phflag=2;

      case 'SU2spin'    % SU(2) total spin

          q=init_qstruct('total spin SU(2)','SU',2);
          q.Sp=sum(SS(:,1));
          q.Sz=sum(SS(:,3));
          if ~isempty(SOP), SOP(end+1)=q; else SOP=q; end

      case 'SUNchannel' % SU(NC) channel symmetry

          if NC<2, error('Wb:ERR',...
            '\n   ERR channel symmetry requires NC>=2 (%g)',NC); end

          Sp=SymOp(1,NC-1); Sz=SymOp(1,NC-1);
          for i=1:NC-1
           % NB! sum over spin!
             Sp(i)=FF(i,1)'*FF(i+1,1) + FF(i,2)'*FF(i+1,2);
           % Sz{i}=0.5*comm(Sp(i),Sp(i)');
             Sz(i)=sum(CZ(1:i))-i*CZ(i+1);
          end
        % nrm=norm(Sp{1},'fro');
        % for i=1:numel(Sz), Sz{i}=Sz{i}*(nrm/norm(Sz{i},'fro')); end
        % keyboard

          q=init_qstruct(sprintf('channel SU(%g)',NC),'SU',NC);
          q.Sp=Sp;
          q.Sz=Sz;

          if ~isempty(SOP), SOP(end+1)=q; else SOP=q; end
          chflag=1; clear Sp Sz

      case 'SpNchannel' % Sp(2*NC) channel symmetry

          if NC<2, error('Wb:ERR',...
             '\n   ERR Sp(2*NC) symmetry requires NC>=2 (%g)',NC); end
           % NB! NC==1 is equivalent to SU2charge!
 
        % build defining representation in Sp(2*NC) for 2*NC-dim spinor psi
        % inheriting all (N-1) Sz and Sp generators from SU(N)
          psi=FF; for i=1:size(FF,1), psi(i,2)=psi(i,2)'; end
          psi=psi(:);

          sz=SymOp(1,NC); sp=SymOp(1,NC);

        % from original defining representation
          for i=1:NC
             is=sprintf('Sp%g:z%g',2*NC,i);
             z=ones(1,NC); if i<NC, z(i+1)=-i; z(i+2:end)=0; end
             sz(i)=SymOp(is,spdiags([z,-z]',0,n2,n2));

             is=sprintf('sp%g:p%g',2*NC,i);
             if i<NC, s=sparse(i,i+1,1,NC,NC);
                  sp(i)=SymOp(is,blkdiag(s,-s'));
             else sp(i)=SymOp(is,sparse(i,i+NC,1,n2,n2)); end
          end

        % transform to irrep order
          if 1
             e=ones(1,NC); e(2:2:end)=-1;
             J=blkdiag(speye(NC),fliplr(spdiags(e',0,NC,NC)));
             sz=transform(sz,J);
             sp=transform(sp,J);

           % NB! must also transform spinor (otherwise *this does
           % not commute with full spin SU(2) symmetry anymore!
           % eg. since Z1 no longer corresponds to proper particl
           % number in some set of channels!) -- Wb,Dec03,11
             psi=psi2psi(J,psi);
           % keyboard
          end

          q=init_qstruct(sprintf('channel Sp(2*%g)',NC),'Sp',2*NC);
          q.Sp=mat2ops(sp,psi);
          q.Sz=mat2ops(sz,psi);

          if ~isempty(SOP), SOP(end+1)=q; else SOP=q; end

          phflag=12; chflag=2;

      otherwise helpthis
      error('Wb:ERR','\n   ERR invalid symmetry %s (%s)',s,Sym_);
    end
  end

  check_commrels(SOP);
  sym=strhcat({SOP.type},'-s',',');

% if 1  % => SU(3): state #1 in IREP-decomp is not MW
        % => problems with inner multiplicity!
% if 0  % better for SU(3)
  if phflag % recommended
   % pre-sort states using diag(Z) [not crucially necessary, though]
   % Wb,Jul11,11
     QZ=catdiag2(cat(2,SOP.Sz));
     [QZ,is]=sortrows(2*QZ); is=flipud(is); QZ=flipud(QZ);

     P=sparse(is,1:D,ones(size(is)), D,D);
     SOP=transform_all(SOP,P); transform(FF,NN,SS,CC,CP,CZ,P);
  else P=[]; end

% get proper symmetry eigenstates
  [U,Is]=getSymStates(strip4mex(SOP)); U=sparse(U);

  SOP=transform_all(SOP,U); transform(FF,NN,SS,CC,CP,CZ,U);
  SOP=resparse(SOP,1E-12);

  if ~isempty(P), U=P*U; end
  
  if nargout>1
     Iout=add2struct('-',NC,SOP,'sym=Sym_',U,Is);
  end

% -----------------------------------------------------------
% F-operators
% -----------------------------------------------------------

  X=FF;
  if phflag, X=[FF,FF']; end % include hermitian conjugate

  k=0; F=QSpace;
  while ~isempty(X), k=k+1;
     [x,Io,X]=getSymmetryOps(X(1),SOP,X);

     for i=1:numel(x), x{i}=full(x{i}); end
     F(k)=compactQS(oc{:},sym,Is.QZ,Is.QZ,Io.QZ, cat(3,x{:}));

  % safeguard: check fermionic commutator relations
     a=QSpace(contractQS(F(k),[1 3],F(k),[1 3],'conjA'));
     b=QSpace(contractQS(F(k),[2 3],F(k),[2 3],'conjB'));
     q=(1/numel(x)) * (a+b);
     if ~isIdentityQS(q), error('Wb:ERR',...
     '\n   ERR fermionic creation/annihilation ops!'); end
  end

  if vflag, % inl(1);
     wblog('<i>','%s',istr); wblog('==>','%s [%s]',Sym_,sym);
     [s,i]=dsize(F(1),'-d');
     s=[ prod(i.sd(:,:),2), permute(prod(i.sc,2),[1 3 2])];
     s=[ sum(s(:)), sum(prod(s,2)) ];
     n=numel(F); if vflag>1, inl(1); end
     wblog(' * ','%g F-op%s @ z=%.3g',n,iff(n~=1,'s',''),s(1)/s(2));
   % if vflag<2, info(F(1)), else display(F(1)); end
     if vflag>1, info(F(1)); end
  end
% keyboard

  F=setOpFlag(F);

% -----------------------------------------------------------
% EZ-operators
% -----------------------------------------------------------

  x=sum(NN(:));
  if norm(x,'-offdiag')>1E-12
     error('Wb:ERR','\n   ERR N-operator not diagonal'); end
  x=full(diag(x.op));
  if norm(x-round(x))>1E-12
     error('Wb:ERR','\n   ERR N-operator is non-integer'); end
  x=round(x); z=zeros(1,size(Is.QZ,2)); oo={};
  E=QSpace(compactQS(oc{:},sym, Is.QZ,Is.QZ,z, eye(D))); oo{end+1}='E';
  Z=QSpace(compactQS(oc{:},sym, Is.QZ,Is.QZ,z, full(diag((-1).^x)))); oo{end+1}='Z';

  if nargout>1
     Iout.E=E; % complete identity space
     Iout.Z=Z; % fermionic signs ( (-1).^(total number of particles) )
  end

% -----------------------------------------------------------
% spin-operators (for every channel and total)
% -----------------------------------------------------------

  S=QSpace;
  if chflag || phflag>10 % SpN symmetry => using Stot only
     SS=sum(SS,1);
  end

  nc=NC/size(SS,1); % Wb,Feb20,13
  if norm(nc-round(nc))>1E-12, error('Wb:ERR',['\n   ERR failed to ' ... 
    'determine effect number of channels for splitup spin operator']);
  end

  for i=1:size(SS,1) % see Stot below
     q=get_spin_ops(SS(i,:),SOP,Is,nc,sym);
     if i>1, S(i,:)=q; else S=q; end
  end

  if vflag
     n=numel(S); if vflag>1, inl(1); end
     wblog(' * ','%g S-op%s',n,iff(n~=1,'s',''));
   % if vflag<2, info(S(1)), else display(S(1)); end
     if vflag>1
       wblog(' * ','%s => %s\N',Sym_,sym); info(S(1));
     end
  end

  if numel(S,1)>1
       Iout.S3=S; S=sum(S,1); % Stot
  else Iout.S3=[]; end

% -----------------------------------------------------------
% operators for particle-hole symmetry
% -----------------------------------------------------------

  if phflag
   % Q = total charge (relative to half-filling)
   % NB! H=U*Q2 (the square is invariant, not Q itself!)

   % every Cz operator is part of an S=1 multiplet!
   % => what is conserved is S.S = (3/4) * (\hat{N}-1)
     n=numel(CZ); if phflag>10, n=1; end
     C2=QSpace(1,n); C3=QSpace(1,n); Q=QSpace;
     for i=1:n
        if phflag>10
           % NB! sum(CZ) is crucial here! CZ(i), while each *fully*
           % exploring the same multiplet, result in a slightly 
           % different basis! // Wb,Dec06,11
             [x,Io]=getSymmetryOps(sum(CZ),SOP);
        else [x,Io]=getSymmetryOps(CZ(i),SOP);
        end
        x=mfull(x); x=cat(3,x{:}); x(find(abs(x)<1E-12))=0;

        C3(i)=compactQS(oc{:},sym, Is.QZ,Is.QZ, Io.QZ, x);
        C2(i)=contractQS(C3(i),[1 3],C3(i),[1 3],'conjA');
     end

     Iout.C3=C3; % spin-op in charge SU(2)
     Iout.C2=C2; % total spin^2 in charge SU(2) space!

     oo{end+1}=sprintf('%gxC3',n);
     oo{end+1}=sprintf('%gxC2',n);

   % total charge
     for i=1:n
     for j=1:n, Q=Q+contractQS(C3(i),[1 3],C3(j),[1 3],'conjA');
     end; end

     Iout.Q2=(4/3)*Q; % total charge
     % factor 4   from (1/2)^2 in definition of Cz
     % factor 1/3 since there are three ops in the IROP C3

     oo{end+1}='Q2';
  end

  if vflag
     wblog(' * ','got { %s } ops',strhcat(oo,'sep',', '));
     if vflag>1, wblog(' * ','%s => %s\N',Sym_,sym); info(Z); end
  end

% -----------------------------------------------------------
% C'C operators (Schrieffer-Wolff) // Wb,Oct03,12
% -----------------------------------------------------------

if Yops
  k=0; Y=QSpace; i=find(CC);

% pre-sort operator set using their z-labels
  [X,Ix]=getSymmetryOps(CC(i),SOP,'-sort');

  try
    while ~isempty(X), k=k+1;
       if vflag, o={'-v'}; else o={}; end
       [x,Io,X]=getSymmetryOps(X(1),SOP,X,o{:},'-P');

       for i=1:numel(x), x{i}=full(x{i}); end
       Y(k)=compactQS(oc{:},sym,Is.QZ,Is.QZ,Io.QZ, cat(3,x{:}));
    end
  catch
    wblog('ERR','failed to obtain Y-ops (SW)'); % Schrieffer-Wolff 
    l=dbstatus; if ~isempty(findcell({l.cond},'error'))
      l=lasterror; disp(l.message); dispstack(l);
      keyboard
    end
  end

  n=numel(Y);
  wblog(' * ','%g Y-op%s (SW)',n,iff(n~=1,'s',''));

  Iout.Y=fixMixedScalar(Y);
end

% wblog('SUC',':)');
% if vflag, inl(1), end
% keyboard

end

% -------------------------------------------------------------------- %
% -------------------------------------------------------------------- %

function S=get_spin_ops(X,SOP,Is,NC,sym)

  e=norm(comm(X(1),X(2))+X(3)); if e>1E-12
  error('Wb:ERR','\n   ERR invalid (sparse) S-ops (%.3g)',e); end

  S=QSpace; k=0; X=fliplr(X); % start getSymmetryOps from z-operators (!)

  while ~isempty(X), k=k+1;
     [x,Io,X]=getSymmetryOps(X(1),SOP,X);

     for j=1:numel(x), x{j}=full(x{j}); end
     S(k)=compactQS(sym,Is.QZ,Is.QZ,Io.QZ, cat(3,x{:}));
  end

% safeguard: check spin commutator relations
% a=QSpace(contractQS(S(1),[2 3],S(2),[1 3],'conjB'));
% q=(1/numel(x)) * (a+b);
% if ~isIdentityQS(q), error('Wb:ERR',...
% '\n   ERR fermionic creation/annihilation ops!'); end

% if numel(S(1).Q)==2 % eg. just got scalar z-operator
  s=zeros(0,3); for i=1:numel(S)
     q=getDimQS(S(1)); q=q(end,:); q(end+1:3)=1;
     s(i,:)=q; % include CGC dimensions!
  end
  if norm(diff(s(:,3)))>1E-12, error('Wb:ERR','invalid S-ops'); end
  if s(1,3)==1 % eg. irreps contain SINGLE operators
     i=reshape(flipud(reshape(1:numel(S),3,[])),size(S));
     S=S(i); % restore original order (Sz comes last)

   % check spin commutator relations
   % NB! must contract S(1) with S(1)', ...

     if numel(S(1).Q)>2
        e=normQS(QSpace(...
          contractQS(S(2),[1 3],S(2),[1 3],'conjA')) ...
        - contractQS(S(1),[1 3],S(1),[1 3],'conjA') - S(3));
     else
        e=normQS(QSpace(...
          contractQS(S(2),1,S(2),1,'conjA')) ...
        - contractQS(S(1),1,S(1),1,'conjA') - S(3));
     end

     if e>1E-12, error('Wb:ERR','\n   ERR invalid spin-CR'); end

   % check total spin
     q=QSpace; for i=1:3
       if numel(S(i).Q)>2    %   V: hermitian conjugate!
            q=q+contractQS(S(i),[1 3],S(i),[1 3],'conjA');
       else q=q+contractQS(S(i),[1  ],S(i),[1  ],'conjA'); end
     end
     e=eigQS(q);
     if norm(max(e(:,1))-(NC/2)*(NC/2+1))>1E-12
     error('Wb:ERR','\n   ERR invalid S-ops'); end

  end

end

% -------------------------------------------------------------------- %
% -------------------------------------------------------------------- %

function [F,Z,Iout]=getLocalSpace_SpinlessFermions(Sym_,varargin)
% function [F,Z,Io]=getLocalSpace_SpinlessFermions(Sym,varargin)
% supported symmetries
%
%   'Acharge'      abelian total charge
%   'SUNchannel'   SU(N) channel symmetry
%
% NB! NO SU(2) particle-hole symmetry in the absence of spin!
% consequently, also no Sp(2n) symmetry
%
% Example: 'Acharge,SUNchannel'
% Wb,Jun27,12

  if nargin<1 || ~ischar(Sym_)
     helpthis, if nargin || nargout
     error('Wb:ERR','invalid usage'), end, return
  end

  getopt('INIT',varargin);
     NC=getopt('NC',1);

     if getopt('-v'), vflag=1; oc={'-v'};
     elseif getopt('-V'), vflag=2; oc={'-V'};
     else vflag=0; oc={}; end

  getopt('check_error');

  istr=sprintf('spinless fermionic %g-level system',NC);

  f2=sparse([0 1; 0 0]); % => [1,0] = empty state
  z2=f2*f2'-f2'*f2; n2=f2'*f2; e2=speye(2);

% build state space
% F-operators
  FF=SymOp(NC,1); NN=SymOp(NC,1);

  for k=1:NC
     q=[ repmat({z2},1,k-1), {f2}, repmat({e2},1,NC-k) ];
     FF(k)=SymOp(sprintf('F%g',k),mkron(q{:}));
     NN(k)=FF(k)'*FF(k); NN(k).istr=sprintf('N%g',k);
  end

  SOP=[]; D=dim(FF);

  Sym=strread(Sym_,'%s','whitespace',';, ')'; clear q

  for s=Sym, s=s{1};
  % change SU<N>channel to generall 'SUNchannel' while checking N==NC
    n=str2num(regexprep(s,'SU(\d+)channel','$1'));
    if ~isempty(n)
       if ~isequal(n,NC), error('Wb:ERR',...
         '\n   ERR invalid symmetry %s (having NC=%g)',s,NC); 
       end
       s='SUNchannel';
    end

    if isequal(s(end-2:end),'(:)') % apply symmetry on each individual channel
         cflag=1; if NC==1, s=s(1:end-3); end % irrelevant
    else cflag=0; end

    switch s
      case 'Acharge'    % abelian total charge
          q=init_qstruct('total charge U(1)','A');
          q.Sz=sum(NN);
          if ~isempty(SOP), SOP(end+1)=q; else SOP=q; end

      case 'Acharge(:)'  % abelian charge for each channel
          for i=1:NC
             q=init_qstruct(sprintf('charge U(1) [%d/%d]',i,NC),'A');
             q.Sz=NN(i);
             if ~isempty(SOP), SOP(end+1)=q; else SOP=q; end
          end

      case 'SUNchannel' % SU(NC) channel symmetry

          if NC<2, error('Wb:ERR',...
            '\n   ERR channel symmetry requires NC>=2 (%g)',NC); end

          Sp=SymOp(1,NC-1); Sz=SymOp(1,NC-1);
          for i=1:NC-1
             Sp(i)=FF(i)'*FF(i+1);
           % Sz{i}=0.5*comm(Sp(i),Sp(i)');
             Sz(i)=0.5*(sum(NN(1:i))-i*NN(i+1));
          end

          q=init_qstruct(sprintf('channel SU(%g)',NC),'SU',NC);
          q.Sp=Sp;
          q.Sz=Sz;

          if ~isempty(SOP), SOP(end+1)=q; else SOP=q; end
          clear Sp Sz

      otherwise helpthis
      error('Wb:ERR','\n   ERR invalid symmetry %s (%s)',s,Sym_);
    end
  end

  check_commrels(SOP);
  sym=strhcat({SOP.type},'-s',',');

  if 1 % recommended
   % pre-sort states using diag(Z)
     QZ=catdiag2(cat(2,SOP.Sz));
     [QZ,is]=sortrows(2*QZ); is=flipud(is); QZ=flipud(QZ);

     P=sparse(is,1:D,ones(size(is)), D,D);
     SOP=transform_all(SOP,P); transform(FF,NN,P);
  else P=[]; end

% get proper symmetry eigenstates
  [U,Is]=getSymStates(strip4mex(SOP)); U=sparse(U);

  SOP=transform_all(SOP,U); transform(FF,NN,U);
  SOP=resparse(SOP,1E-12);

  if ~isempty(P), U=P*U; end

  if nargout>1, Iout=add2struct('-',NC,SOP,'sym=Sym_',U,Is); end

% -----------------------------------------------------------
% F-operators
% -----------------------------------------------------------

  X=FF;

  k=0; F=QSpace;
  while ~isempty(X), k=k+1;
     [x,Io,X]=getSymmetryOps(X(1),SOP,X);

     for i=1:numel(x), x{i}=full(x{i}); end
     F(k)=compactQS(oc{:},sym,Is.QZ,Is.QZ,Io.QZ, cat(3,x{:}));

  % safeguard: check fermionic commutator relations
     a=QSpace(contractQS(F(k),[1 3],F(k),[1 3],'conjA'));
     b=QSpace(contractQS(F(k),[2 3],F(k),[2 3],'conjB'));
     q=(1/numel(x)) * (a+b);
     if ~isIdentityQS(q), error('Wb:ERR',...
     '\n   ERR fermionic creation/annihilation ops!'); end
  end

  if vflag, % inl(1);
     wblog('<i>','%s',istr); wblog('==>','%s [%s]',Sym_,sym);
     [s,i]=dsize(F(1),'-d');
     s=[ prod(i.sd(:,:),2), permute(prod(i.sc,2),[1 3 2])];
     s=[ sum(s(:)), sum(prod(s,2)) ];
     n=numel(F); if vflag>1, inl(1); end
     wblog(' * ','%g F-op%s @ z=%.3g',n,iff(n~=1,'s',''),s(1)/s(2));
   % if vflag<2, info(F(1)), else display(F(1)); end
     if vflag>1, info(F(1)); end
  end
% keyboard

  F=setOpFlag(F);

% -----------------------------------------------------------
% E,Q,Z-operators
% -----------------------------------------------------------

  x=sum(NN(:));
  if norm(x,'-offdiag')>1E-12
     error('Wb:ERR','\n   ERR N-operator not diagonal'); end
  x=full(diag(x.op));
  if norm(x-round(x))>1E-12
     error('Wb:ERR','\n   ERR N-operator is non-integer'); end
  x=round(x); z=zeros(1,size(Is.QZ,2)); oo={};
  E=QSpace(compactQS(oc{:},sym, Is.QZ,Is.QZ,z, eye(D))); oo{end+1}='E';
  Z=QSpace(compactQS(oc{:},sym, Is.QZ,Is.QZ,z, full(diag((-1).^x)))); oo{end+1}='Z';

  if nargout>1
     Iout.E=E; % complete identity space
     Iout.Z=Z; % fermionic signs ( (-1).^(total number of particles) )
  end

  if vflag
     wblog(' * ','got { %s } ops',strhcat(oo,'sep',', '));
     if vflag>1, wblog(' * ','%s => %s\N',Sym_,sym); info(Z); end
  end

end

% -------------------------------------------------------------------- %
% -------------------------------------------------------------------- %

function [S,Iout]=getLocalSpace_Spin(qloc,varargin)
% function [S,Io]=getLocalSpace_Spin(qloc,varargin)
% Wb,Jul30,12

  if nargin<1 || ~isnumber(qloc)
     helpthis, if nargin || nargout
     error('Wb:ERR','invalid usage'), end, return
  end

  getopt('INIT',varargin);
     Aflag=getopt('-A'); % use abelian Sz only

     if getopt('-v'), vflag=1; oc={'-v'};
     elseif getopt('-V'), vflag=2; oc={'-V'};
     else vflag=0; oc={}; end

  getopt('check_error');

  istr=sprintf([iff(Aflag,'abelian ',''), 'spin-%g system'],qloc);

  if norm(round(2*qloc)-2*qloc) || qloc<=0
  error('Wb:ERR','\n   ERR invalid spin S=%g',qloc); end

% NB! use -sym flag, such that also in abelian case S^2 = S.S (!)
  sx=spinmat(2*qloc+1,'-sp','-sym');
  sx=sx([1 3 2]); % [ Sp, Sm, Sz ] order

  SS=SymOp(1,3); clear q
  for j=1:3
     s=getOpName('spin',1,j);
     SS(1,j)=SymOp(s([1,3:end]),sx{j});
  end

  if Aflag
     q=init_qstruct('total spin U(1)','A');
     q.Sz=2*sum(SS(:,3));
  else
     q=init_qstruct('total spin SU(2)','SU',2);
     q.Sp=sum(SS(:,1));
     q.Sz=sum(SS(:,3));
  end
  SOP=q; D=dim(SS);

  check_commrels(SOP);
  sym=strhcat({SOP.type},'-s',',');

  if 1 % recommended
   % pre-sort states using diag(Z)
     QZ=catdiag2(cat(2,SOP.Sz));
     [QZ,is]=sortrows(2*QZ); is=flipud(is); QZ=flipud(QZ);

     P=sparse(is,1:D,ones(size(is)), D,D);
     SOP=transform_all(SOP,P); transform(SS,P);
  else P=[]; end

% get proper symmetry eigenstates
  [U,Is]=getSymStates(strip4mex(SOP)); U=sparse(U);

  SOP=transform_all(SOP,U); transform(SS,U);
  SOP=resparse(SOP,1E-12);
  
  if ~isempty(P), U=P*U; end

  if nargout>1
     Iout=add2struct('-','Sloc=qloc',SOP,'sym=''SpinS''',U,Is);
  end

% -----------------------------------------------------------
% further operators: E
% -----------------------------------------------------------

  z=zeros(1,size(Is.QZ,2)); oo={};
  E=QSpace(compactQS(oc{:},sym, Is.QZ,Is.QZ,z, eye(D))); oo{end+1}='E';

  if nargout>1
     Iout.E=E; % complete identity space
  end

  if vflag
     wblog(' * ','got { %s } ops',strhcat(oo,'sep',', '));
     if vflag>1, wblog(' * ','%s => %s\N',Sym_,sym); end
  end

% -----------------------------------------------------------
% spin-operators
% -----------------------------------------------------------

  k=0; S=QSpace; X=sum(SS,1); oo={};

  e=norm(comm(X(1),X(2))+X(3)); if e>1E-12
  error('Wb:ERR','\n   ERR invalid (sparse) S-ops (%.3g)',e); end

  X=fliplr(X); % start getSymmetryOps from z-operators (!)

  while ~isempty(X), k=k+1;
     [x,Io,X]=getSymmetryOps(X(1),SOP,X);

     for i=1:numel(x), x{i}=full(x{i}); end
     S(k)=compactQS(oc{:},sym,Is.QZ,Is.QZ,Io.QZ, cat(3,x{:}));
  end

  S=fixMixedScalar(S);

  if qloc==1, oo{end+1}='and T-op';
   % get T-operator: this is the 2nd non-trivial other IROP having S=2,
   % other than S3 having S=1; follow up to discussion with Wei Li.
   % Wb,Nov03,12
     [x,Io]=getSymmetryOps(SS(1)*SS(1),SOP);
     for i=1:numel(x), x{i}=full(x{i}); end
     Iout.T=compactQS(oc{:},sym,Is.QZ,Is.QZ,Io.QZ, cat(3,x{:}));
  end

% check total spin
  if numel(S)==1
     q=contractQS(S,[1 3],S,[1 3],'conjA');
  else
     X=QSpace;
     for i=1:numel(S)
        if numel(S(i).Q)>2
             q=contractQS(S(i),[1 3],S(i),[1 3],'conjA');
        else q=contractQS(S(i),[1  ],S(i),[1  ],'conjA'); end
        X=X+q;
     end
     q=X;
  end

  e=eigQS(q);
  if norm(max(e(:,1))-qloc*(qloc+1))>1E-12
  error('Wb:ERR','\n   ERR invalid S-ops'); end

  if vflag
     Iout.S2=QSpace(q);
     n=numel(S); if vflag>1, inl(1); end
     wblog(' * ',['%g S-op%s ' oo{:}],n,iff(n~=1,'s',''));
   % if vflag<2, info(S(1)), else display(S(1)); end
     if vflag>1
       wblog(' * ','%s => %s\N',Sym_,sym); info(S(1));
     end
  end

end

% -------------------------------------------------------------------- %
% -------------------------------------------------------------------- %

function [S,Iout]=getLocalSpace_SpinSUN(sym,qloc,varargin)
% function [S,Iout]=getLocalSpace_SpinSUN(sym,qloc [,opts])
%
%    qloc is the symmetry of the local multiplet to choose
%    (default: defining representation).
%
%    NB! qloc needs to be specified in terms of general SU(N)
%    notation based on the underlying Young diagram, hence
%    qloc=2*S in the case of SU(2) [note that this is different
%    from getLocalSpace_Spin which interprets qloc in terms
%    of standard half-spins].
%
% Wb,Apr07,14

   if nargin<2, helpthis
      error('Wb:ERR','\n   ERR invalid usage'); end
   if isempty(regexp(sym,'^SU\d+$')), helpthis
      error('Wb:ERR','\n   ERR invalid symmetry ''%s''',sym); end
   if ~isnumeric(qloc), helpthis, qloc
      error('Wb:ERR','\n   ERR invalid qloc for symmetry ''%s''',sym);
   end

   getopt('INIT',varargin);
      if getopt('-v'), vflag=1; oc={'-v'};
      elseif getopt('-V'), vflag=2; oc={'-V'};
      else vflag=0; oc={}; end
   getopt('check_error');

   N=str2num(sym(3:end)); r=N-1;
   [F,Z,I]=getLocalSpace('Fermion','SUNchannel','NC',N,oc{:});

   I.SOP(1).info=regexprep(I.SOP(1).info,'channel','Spin');

   if isempty(qloc), qloc=[1, zeros(1,r-1)];
   elseif numel(qloc)~=r || size(qloc,1)~=1 || ~isequal(round(qloc),qloc)
      qloc, error('Wb:ERR','\n   ERR invalid qloc');
   end

 % define spin operator through "Schrieffer-Wolff transformation"
   q1=F.Q{3}(1,:); % dual to defining representation ("annihilation op")
   i=matchIndex(I.E.Q{1},q1);
   if isempty(i), error('Wb:ERR',...
      '\n   ERR failed to identify defining representation !?');
   end

 % e.g. qloc=[020] // Thomas Quella
   A0=QSpace(getIdentityQS(I.E,1,I.E,1));
   i=matchIndex([A0.Q{1} A0.Q{2}],[fliplr(q1), q1]);
   q=A0.Q{3}(i,:); j=find(sum(q~=0,2));
   if numel(i)~=2 || numel(j)~=1, error('Wb:ERR',...
      '\n   ERR failed to find regular representation ("spin")');
   end
   qs=q(j,:); % multiplet label of "spin operator" (reg. representation)
   is=i(j);

   istr=sprintf('%s-spin operator ([%s])',sym,vec2str(qs));

   i=matchIndex(I.E.Q{1}, qloc);
   if ~isempty(i) % Wb,Feb25,15
      Eq=getsub(I.E,i);
      A0=QSpace(getIdentityQS(Eq,1,Eq,1));
      Ak=A0;
      for k=0:4
         i=matchIndex(Ak.Q{3}, qs); if ~isempty(i), break; end
         Ak=QSpace(getIdentityQS(Ak,3,Eq,1));
      end
      if isempty(i), error('Wb:ERR',...
        '\n   ERR failed to generate spin multiplet'); end
      Ak=getsub(Ak,i);
      Ak=QSpace(getIdentityQS(Ak,2,Ak,3));

      i=matchIndex(Ak.Q{3}, qloc);
      if numel(i)~=1, error('Wb:ERR',...
        '\n   ERR failed to generate spin operator'); end

    % NB! require qdir=+-- for operators (!) // Wb,Jan14,15
      S=QSpace(permuteQS(getsub(Ak,i),[3 1 2],'conj'));

    % normalize S2 to one for now // fix normalization later
      S.data{1}=1;
      S2=QSpace(contractQS(S,'13*',S,'13'));
      S.data{1}=1/sqrt(S2.data{1}(1));

      d=double(S.info.cgr(1).size);
    % keyboard, return
   else
    % NB! qs must be part of "input" spaces (otherwise "arrows"
    % are messed up!). The spin (=regular) representation, however,
    % is NEVER part of I.E! => need a 2nd round // Wb,Apr23,14

      if ~isequal(qloc,fliplr(qloc)), qloc=[qloc; fliplr(qloc)]; end
      [i1,i2,Im]=matchIndex(A0.Q{3},qloc);
      if isempty(i1), error('Wb:ERR',...
         '\n   ERR failed to find qloc in state space decomposition');
      end

      A1=QSpace(getIdentityQS(getsub(A0,i1),3,getsub(A0,is(1)),3));
      i=matchIndex([A1.Q{1}, A1.Q{2}, A1.Q{3}],...
        [qloc(end,:), qs, qloc(1,:)]);
    % S=permute(getsub(A1,i),[1 3 2]);
    % NB! require qdir=+-- for operators (!) // Wb,Jan14,15
      S=QSpace(permuteQS(getsub(A1,i),[3 1 2],'conj'));

      if numel(i)>1
         error('Wb:ERR',['\n   ' ... 
           'ERR failed to identify spin operator (%s [%s], qs=[%s])\n   ' ...
           'ERR due to outer multiplicity [%s] * [%s] => [%s](%g)\n'], ...
         sym, vec2str(qs), vec2str(qloc(1,:)), ...
         vec2str(qloc), vec2str(qs), vec2str(qloc), numel(i));
      elseif isempty(i)
         if isequal(sym,'SU2'), h=['\n   ERR ==> hint: ' ... 
           'having SU2, use getLocalSpace(''Spin'',S,...), instead.\n'];
         else h=''; end

         error('Wb:ERR',...
         ['\n   ERR failed to find spin representation (%s [%s], qs=[%s])' h],...
         sym, vec2str(qs), vec2str(qloc(1,:)));
      end

    % choose normalization of spin operator consistent with SU(2)
    % (the generalization to arbitrary N is a bit heuristic, though,
    % yet note that it has to be a quadratic polynomial in qloc).
    % Wb,Apr07,14
      d=double(S.info.cgr(1).size);
      S.data{1}=sqrt(d(2)*(qloc*(qloc+2)')/4) / ...
         norm(mpfr2dec(S.info.cgr(1).cgw),'fro');
       % NB! relying on CStore/CData convention that |cdata|^2=1
   end

   S.info.otype='operator';

   if nargout<2, return; end

   Iout.istr=istr;
   q=add2struct(rmfield(I,{'Z','NC'}),qloc,qs,A0,'A1?');

   Iout=structmerge(Iout,q);
   Iout.sym=sprintf('Spin%s',sym);

   if vflag
      wblog('<i>','setup for %s (d=%d)',istr,d(end));
      wblog(' * ','having qloc=[%s] (d=%d)', vec2str(qloc),d(1));
   end

 % further operators: E

 % Iout.E=QSpace(contractQS(S,[1 3],S,[1 3],'conjA'));
 % Iout.E.data{1}=1;

 % keyboard

end

% -------------------------------------------------------------------- %
% -------------------------------------------------------------------- %

function [S,Iout]=getLocalSpace_SUN(N,varargin)
% function [S,Io]=getLocalSpace_SUN(N,varargin)
%    plain SU(N) site with single local IROP,
%    corresponding to the regular represenation.
% Wb,Nov03,12

  if nargin<1 || ~isnumber(N)
     helpthis, if nargin || nargout
     error('Wb:ERR','invalid usage'), end, return
  end

  getopt('INIT',varargin);
     if getopt('-v'), vflag=1; oc={'-v'};
     elseif getopt('-V'), vflag=2; oc={'-V'};
     else vflag=0; oc={}; end
  getopt('check_error');

  istr=sprintf('SU(%g) site',N);

  if N<2 || N>8, error('Wb:ERR',...
    '\n   ERR SU(%g) not yet implemented !?',N);
  end

  Sp=SymOp(1,N-1); Sz=SymOp(1,N-1);
  for i=1:N-1
     Sp(i)=SymOp(sprintf('Sp%g',i),sparse(i,i+1,1,N,N));
     z=ones(N,1); z(i+1)=-i; z(i+2:end)=0;
     Sz(i)=SymOp(sprintf('Sz%g',i),spdiags(z/2,0,N,N));
  end

  q=init_qstruct(sprintf('plain SU(%g) site',N),'SU',N);
  q.Sp=Sp;
  q.Sz=Sz;

  SOP=q; D=dim(Sp); clear Sp Sz
  check_commrels(SOP);
  sym=strhcat({SOP.type},'-s',',');

  if 1 % recommended
   % pre-sort states using diag(Z)
     QZ=catdiag2(cat(2,SOP.Sz));
     [QZ,is]=sortrows(2*QZ); is=flipud(is); QZ=flipud(QZ);

     P=sparse(is,1:D,ones(size(is)), D,D);
     SOP=transform_all(SOP,P); % transform(SS,P);
  else P=[]; end

% get proper symmetry eigenstates
  [U,Is]=getSymStates(strip4mex(SOP)); U=sparse(U);

  SOP=transform_all(SOP,U); % transform(SS,U);
  SOP=resparse(SOP,1E-12);
  
  if ~isempty(P), U=P*U; end

  if nargout>1
     Iout=add2struct('-',N,SOP,sprintf('sym=''SU%g''',N),U,Is);
  end

  wblog('<i>','setting up %s',istr);

% -----------------------------------------------------------
% get operators
% -----------------------------------------------------------

% identity operator E
  z=zeros(1,size(Is.QZ,2)); oo={};
  E=QSpace(compactQS(oc{:},sym, Is.QZ,Is.QZ,z, eye(D))); oo{end+1}='E';

  if nargout>1
     Iout.E=E; % complete identity space
  end

  if vflag
     s=iff(numel(oo)>1,'s','');
     wblog(' * ','got { %s } op%s',strhcat(oo,'sep',', '),s);
  end

% get IROP for regular representation
% (adjoint representation = spin operator)
  [x,Io]=getSymmetryOps(SymOp('MW',sparse(1,N,1,N,N)),SOP);
  for i=1:numel(x), x{i}=full(x{i}); end
  S=QSpace(compactQS(oc{:},sym,Is.QZ,Is.QZ,Io.QZ, cat(3,x{:})));

  if vflag
     wblog(' * ','got IROP for regular representation');
  end

end

% -------------------------------------------------------------------- %
% -------------------------------------------------------------------- %
% keep scalar operator dimension (3rd index)
% if non-scalar operators exist within the same QSpace array

function F=fixMixedScalar(F)

   if numel(F)<2, return; end

   mark=zeros(size(F)); n=numel(F);
   for i=1:n
      if numel(F(i).Q)==3, q=norm(F(i).Q{3});
         if q>1E-12, mark(i)=1; else mark(i)=-1; end
      end
   end

   if any(mark(:)>0)
      I=find(mark==0); n=numel(I);
      for k=1:n, i=I(k);
          F(i).Q{3}=zeros(size(F(i).Q{1}));
      end
   else
      I=find(mark<0); n=numel(I);
      for i=1:n, F(i)=squeeze(F(i)); end
   end

end

% -------------------------------------------------------------------- %
% -------------------------------------------------------------------- %

function s=getOpName(tag,varargin)

   switch tag
      case 'fops'

         if numel(varargin)~=2
            error('Wb:ERR','\n   ERR invalid usage (%s)',tag); end
         [i,s]=deal(varargin{:});

         switch s
            case 1, s='u';
            case 2, s='d';
            otherwise s, error('Wb:ERR','\n   ERR invalid spin');
         end
         s=sprintf('F%g%s',i,s);

      case 'spin'
         if numel(varargin)~=2
            error('Wb:ERR','\n   ERR invalid usage (%s)',tag); end
         [i,s]=deal(varargin{:});

         switch s
            case 1, s='+';
            case 2, s='-';
            case 3, s='z';
            otherwise s, error('Wb:ERR','\n   ERR invalid spin');
         end
         s=sprintf('S%g(%s)',i,s);

      otherwise, tag, error('Wb:ERR','\n   ERR invalid tag'); 
   end
end

% -------------------------------------------------------------------- %
% -------------------------------------------------------------------- %

function F=setOpFlag(F)

   for i=1:numel(F)
      if numel(F(i).Q)==3,F(i).info.otype='operator'; end
   end

end

% -------------------------------------------------------------------- %
% -------------------------------------------------------------------- %

function SOP=transform_all(SOP,U)

   s=inputname(2); if ~isempty(s), o={'-n',s}; else o={}; end
   for i=1:numel(SOP)
      SOP(i).Sp=transform(SOP(i).Sp,o{:},U);
      SOP(i).Sz=transform(SOP(i).Sz,o{:},U);
   end

end

function SOP=resparse(SOP,varargin)

   for i=1:numel(SOP)
      SOP(i).Sp=sparse(SOP(i).Sp, varargin{:});
      SOP(i).Sz=sparse(SOP(i).Sz, varargin{:});
   end

end

% prepare SOP to mex-file-input such as getSymStates
% Wb,Dec05,11
function SOP=strip4mex(SOP)

   for i=1:numel(SOP)
      SOP(i).Sp=getops(SOP(i).Sp,'-f');
      SOP(i).Sz=getops(SOP(i).Sz,'-f');
   end

end

% -------------------------------------------------------------------- %
% -------------------------------------------------------------------- %
% Wb,Dec03,11

function G=mat2ops(G,psi)

  s=size(psi(1)); q0=sparse(s(1),s(2));

  for k=1:numel(G)
     [i,j,g]=find(G(k).op); q=SymOp('',q0);
     for l=1:numel(g), q=q+g(l)*psi(i(l))'*psi(j(l)); end
     G(k)=q;
  end

end

function Upsi=psi2psi(U,psi)
% |psi> => U|psi>

  if numel(psi)~=size(U,1), error('Wb:ERR',...
    '\n   ERR size mismatch (%g/%g)',size(U,1),numel(psi)); end
  e=norm(U*U'-speye(size(U)),'fro'); if e>1E-12, error('Wb:ERR',...
    '\n   ERR unitary matrix U expected! (%.3g)',e); end

  s=dim(psi); Upsi=SymOp(size(psi)); q0=sparse(s,s);

  for i=1:numel(psi)
     [x,j,g]=find(U(i,:)); q=SymOp('',q0);
     for l=1:numel(g), q=q+g(l)*psi(j(l)); end
     Upsi(i)=q;
  end
end

% -------------------------------------------------------------------- %
% -------------------------------------------------------------------- %
% check all commutators for symmetry operations (within sparse still)
% outsourced from getLocalSpace_SpinfullFermions()
% Wb,Jul31,12

function check_commrels(SOP)

  for i=1:numel(SOP)
     S1=SOP(i).Sz; n1=numel(S1);
     for i1=1:n1, for i2=i1+1:n1
        e=norm(comm(S1(i1),S1(i2)));
        if e>1E-12, error('Wb:ERR','\n   ERR invalid z-operator'); end
     end, end

     S1=[ SOP(i).Sz, SOP(i).Sp ]; n1=numel(S1);

     for j=i+1:numel(SOP)
        S2=[ SOP(j).Sz, SOP(j).Sp ]; n2=numel(S2);
        for i1=1:n1, for i2=1:n2
           e=norm(comm(S1(i1),S2(i2)));
           if e>1E-12, inl(1); disp(SOP(i)), disp(SOP(j))
              error('Wb:ERR',['\n   ERR incompatible symmetries ' ...
             'SOP(%d).%d and SOP(%d).%d (e=%g)'],i,i1,j,i2,e);
           end
        end, end
     end
  end
end

% -------------------------------------------------------------------- %
% -------------------------------------------------------------------- %

