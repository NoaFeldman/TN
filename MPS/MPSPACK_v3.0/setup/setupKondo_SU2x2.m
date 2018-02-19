% -------------------------------------------------------------------- %
% setup operator set for standard Kondo problem (single channel)
% using particle-hole and spin SU(2) symmetry (if B=0)
% Wb,Sep19,11
% -------------------------------------------------------------------- %

% adapted from setupSIAM.m

  global param

  if ~exist('J','var')
     J=0.12; B=0; Lambda=2; N=50;
  else
     setdef('B',0);
     setdef('Lambda',2);
  end

  istr=''; sym='';
  param=add2struct('-',istr,sym,'fJ=J',B,Lambda,N);

% pg_r = pseudo-gag coefficient r
  co={'-x'};
  if exist('z',   'var'), co(end+1:end+2)={'z',z   }; param.z=z;    end
  if exist('pg_r','var'), co(end+1:end+2)={'r',pg_r}; param.r=pg_r; end
  setopts(co,'Nmax?');

  if exist('ALambda','var')
     if ischar(ALambda) || isempty(ALambda) || ALambda
        co(end+1)={'-AL'}; param.ALambda=1;
     else
        co(end+1)={'-w'}; % param.ALambda=0;
        if isfield(param,'ALambda') param=rmfield(param,'ALambda'); end
     end
  end

  [ff,q,Ic]=getNRGcoupling(1,Lambda,N,co{:});

  if exist('J','var')
       ff(1)=J;
  else ff(1)=8*Gamma/(pi*U); end

  param.fJ=ff(1);

  if exist('TK0')==1, TK=TK0; else
  TK=TKondo; end

% -------------------------------------------------------------------- %
% operator setup
% -------------------------------------------------------------------- %

  setdef('NC',1);
  param.NC=NC; o={'NC',NC,'-v'};

  param.istr=sprintf('%s-channel Kondo(J) model',num2text(NC));

  if NC>1
     if 1 % Wb,Feb02,14
        if ~isset('B') % Wb,Aug20,13
             sym='SU2spin,SpNchannel';
        else sym='Aspin,SpNchannel'; end
     else
        if ~isset('B') % Wb,Aug20,13
             sym='Acharge,SU2spin,SUNchannel';
        else sym='Acharge,Aspin,SUNchannel'; end
     end

     [FF,Z,SS,IS]=getLocalSpace('FermionS',sym,o{:});

  elseif isset('SYM') && isequal(SYM,'A,A')
     [FF,Z,SS,IS]=getLocalSpace('FermionS','Acharge,Aspin',  o{:});
  elseif isset('SYM') && isequal(SYM,'SU2,A') || isset('B') % Wb,Sep28,13
     [FF,Z,SS,IS]=getLocalSpace('FermionS','SU2charge,Aspin',  o{:});
  elseif isset('SYM') && isequal(SYM,'A,SU2')
     error('Wb:ERR','\n   ERR invalid SYM=%s for Kondo model !??',SYM);
   % NB! Kondo model is always particle-hole symmetric!
   % => no A,SU2 symmetry!
  else
     [FF,Z,SS,IS]=getLocalSpace('FermionS','SU2charge,SU2spin',o{:});
   % S0=SS; for i=1:numel(S0), S0(i)=contractQS(IS.E,2,S0(i),1); end
  end

  SYM=getsym(FF); param.sym=IS.sym;

  E1=QSpace(getIdentityQS(SS(end)));
  % eliminate irrelevant (decoupled) space from impurity!
  % Wb,May09,14

  A0=QSpace(permuteQS(getIdentityQS(E1,1,Z,1),[1 3 2]));
  A1=QSpace(permuteQS(getIdentityQS(A0,2,Z,1),[1 3 2]));

% NB! for consistency with non abelian setup: (2J)S.s = J Sigma.s
  f=2*ff(1);

  if numel(SS)>1, H0=QSpace;
     for i=1:numel(SS), q=SS(i);
        if rank(q)>2
           q=contractQS(...
              contractQS(q,2,A0,1),[1 4 2],...  % LoRs
              contractQS(A0,3,q,2),[1 3 4],...  % LRso
             'conjA' ...
           );
        else
           q=contractQS(...
              contractQS(q,2,A0,1),[1 3],... % LRs
              contractQS(A0,3,q,2),[1 3],... % LRs
             'conjA' ...
           );
        end
        H0=H0+q;
     end
     H0=f*H0;
  else
     H0=f*QSpace(contractQS(...
        contractQS(SS,2,A0,1),[1 4 2],...  % LoRs
        contractQS(A0,3,SS,2),[1 3 4],...  % LRso
       'conjA') ...
     );
  end
  
  ff=ff(2:end); % done with the first iteration!!

% -------------------------------------------------------------------- %
% apply B field // Wb,Jul13,13

  if isset('B')
     if numel(SS)<3 || dim(SS(end),'-op')~=1, error('Wb:ERR',...
       '\n   ERR invalid symmetry/spin setting for B!=0');
     end
     q=SS(end); if rank(q)>2
        if dim(q,'-op')~=1, error('Wb:ERR',...
           '\n   ERR unexpected rank-%g spin operator',rank(q)); end
        q=fixScalarOp('-f',q);
     end
   % Wb,Aug21,13: fixed B contribution (hint M. Hanl)
     q=QSpace(contractQS(contractQS(q,2,A0,1),[1 3],A0,[1 3],'conjA'));
   % NB! use H=-B*Sz to obtain positive <Sz> for B>0!
     H0=H0-B*q;
  end

% -------------------------------------------------------------------- %

% H0=H0+0*QSpace(contractQS(A0,[1 3],A0,[1 3],'conjA'));
% irrelevant by now (plus calls skipzeros!) // Wb,Jul13,13

% NB! NO LONGER NEEDED (NRG::KEEP_AH0) // Wb,Apr13,13
% [a,I]=eigQS(H0);
% Ax=QSpace(contractQS(A0,2,I.AK,1,[1 3 2]));
% Hx=diag(QSpace(I.EK));

  op1=[]; op2=FF;
% spectral function (composite operator for Kondo model)
  if isset('B'), m=numel(FF); else m=1; end
  for i=1:m
     q=QSpace(contractQS(A0,[1 3],contractQS(A0,3,FF(i),2),[1 3],'conjA'));
     op2(i)=QSpace(contractQS(H0,2,q,1)) - contractQS(H0,1,q,2,[2 1 3]);
     op2(i).info.otype=FF(i).info.otype;
  end
  i=m+1;

% FC=FF;
  FC=op2; % Wb,Sep28,13

% spin-spin correlation function for magnetic susceptibility
  q=SS(end); if rank(q)>2, p={[1 3 4 2]}; else p={}; end
  op2(i)=QSpace(contractQS(A0,[1 3],contractQS(q,2,A0,1,p{:}),[1 3],'conjA'));
  op2(i).info.otype=q.info.otype;

  zflags=ones(1,m); zflags(m+1)=0;

  Z0=QSpace(contractQS(A0,[1 3],contractQS(A0,3,Z,2),[1 3],'conjA'));
% keep for general symmetries! (even in presence of particle-hole
% alternation, where application of Z0 with the impurity should
% be fine with ZFLAGS=2! [cf NRGWilsonQS])

% -------------------------------------------------------------------- %

  if SYM(1)=='S', ZFLAG=2; else clear ZFLAG; end
  if isequal(SYM,'A,A')
     fixAbelian(FF,FC,SS,op1,op2,'-xop3');
  end

  wblog('<i>', 'Kondo parameters (TK=%.4e)\N', TK);
  disp(param);

% -------------------------------------------------------------------- %

