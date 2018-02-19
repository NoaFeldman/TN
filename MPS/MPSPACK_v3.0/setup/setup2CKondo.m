
% -------------------------------------------------------------------- %
% operator setup for two-channel Kondo problem
% NB! see also setupKondo with NC=2 (!)
% NB! see also setupKondo_SU2x2.m for 1-channel setup
% Wb,Mar01,08
% Wb,Apr24,12 : adapted to non-abelian symmetries [SU(2)_spin * Sp(4)]
% -------------------------------------------------------------------- %

% NC=2; % two channels (modified from setupKondo_AM)

% NB! J=0.12 is TOO SMALL (TK<1E-8)
  setdef(...
     'J', 0.25, 'B',0E-6, 'NC',2, ...
     'N',36, 'Lambda',4,'Nkeep',1024,'Etrunc',8 ...
  );

  global param
  param=add2struct('-',J,B,Lambda,N,Nkeep,Etrunc);
  TK=TKondo;

  if exist('z','var'), co={'z',z}; param.z=z; else co={}; end
  ff=getNRGcoupling(1,Lambda,N,'-x',co{:});
  ff=ff(2:end);

  wblog('<i>', 'Kondo parameters (NC=%g, TK=%.4g)\N',NC,TK);

% possible charge/channel symmetries // choose by setting SYM

  setdef('SYM',''); % assume spin symmetry comes first
  if B || ~isempty(regexp(SYM,'^A,')); % Bflag
       sym={'Aspin',   'SpNchannel'};
  else sym={'SU2spin', 'SpNchannel'};
  end

  if ~isempty(SYM)
     if     ~isempty(regexp(SYM,',Sp'))       sym{2}='SpNchannel';
     elseif ~isempty(regexp(SYM,',A,SU'))     sym{2}='Acharge,SUNchannel';
     elseif ~isempty(regexp(SYM,',A,A'))      sym{2}='Acharge(:)';
     elseif ~isempty(regexp(SYM,',SU2,SU2'))  sym{2}='SU2charge(:)';
     end
  end

% NB! rnrg.m uses FC
  [~,~, s0, ~]=getLocalSpace('FermionS', sym{1},'NC',1);
  [FC,Z,SS,IS]=getLocalSpace('FermionS',[sym{1} ',' sym{2}],'NC',NC,'-v');

% use same symmetry setting in latter calls (unless cleared)
  SYM=getsym(IS.E); param.sym=IS.sym;

  disp(param);

  if ~isempty(regexp(sym{2},'^(SpN|SU2)'))
     % NB! Sp4 contains particle-hole symmetry!
     % 2nd Wilson site below is added *without* Z operator
     % => local site of A0 has its fermionic order flipped
     % => ZFLAG=3 // ok. Wb,May05,14
       ZFLAG=3;
  else ZFLAG=1; end

  s0_=s0;

% [d,s]=getsym(FC(1),'-d',2);
% s0=appendScalarSymmetry(s0,s); % Wb,Feb06,14
  s=getsym(FC(1),'-c'); % s0_=s0;
  for i=2:numel(s)
     s0=appendScalarSymmetry(s0,s{i}); % Wb,Mar14,14
  end

  A0=QSpace(permuteQS(getIdentityQS(s0(end),Z),[1 3 2]));

  r=numel(s0); % e.g. r=3 if B!=0 (Sx,Sy,Sz)
  if numel(SS)~=r || (r~=1 && r~=3) % Wb,Apr26,16
     error('Wb:ERR','\n   ERR invalid spin setting'); end

  H0=QSpace;
  for i=1:r, Q=contract(J*s0(i),'1*',A0,1);
     if numel(s0(i).Q)==3
          Q=contract(Q,'42',SS(i),'23');
     else Q=contract(Q,3,SS(i),2);
     end
     Q=contract(A0,'13*',Q,'13');
     H0=H0+Q;
  end
  HJ=H0;

  if B % Bflag
   % -B*Sz to ensure postive magnetization <Sz> // Wb,Apr28,16
     H0=H0+contract(A0,'13*',contractQS(-B*s0(end),2,A0,1),'13');
  end

% NB! composite operators for T-matrix derive from commutator [F,HJ]
% see PHYS-notes
  FX=FC; OX=FC;
  for i=1:numel(FC)
     FX(i)=contract(A0,'13*',contractQS(A0,3,FC(i),2),'13');
     OX(i)=QSpace(contractQS(FX(i),2,HJ,1,[1 3 2])) - contractQS(HJ,2,FX(i),1);
  end

% add operator for spin susceptibility
  op2=OX; op1=[];
  op2(end+1)=FX(1);
% op2(end+1)=contractQS(A0,'13*',contractQS(A0,3,SS(end),2),'13');
  op2(end+1)=contractQS(A0,'13*',contractQS(A0,1,s0(end),2),'32');

  for i=1:numel(op2)
     if numel(op2(i).Q)>2
     op2(i).info.otype='operator'; end
  end
  zflags=ones(1,numel(op2)); zflags(end)=0;

  Z0=contract(A0,'13*',contractQS(A0,3,Z,2),'13');

% add 2nd Wilson site for simple treatment of composite operators
  AJ=A0; H0_=H0;
  A0=QSpace(permuteQS(getIdentityQS(Z0,Z),[1 3 2]));

  f1=ff(1); ff=ff(2:end);

  H0=contract(A0,'13*',contractQS(H0,2,A0,1),'13');

  for i=1:numel(FC)
     if ZFLAG>1
      % NO (Z*F) here (NB! Sp4 also contains particle-hole symmetry)
        H0=H0+contract(A0,'13*',...
           contractQS(contractQS(FX(i),'1*',A0,1),'42',...
           f1*FC(i),'23'), ...
        '13');
     else
      % regular treatment of fermionic signs using Z // Wb,Mar14,14
        Q=contract(A0,'13*',...
           contractQS(contractQS(FX(i),'1*',A0,1),'42',...
           contractQS(f1*Z,2,FC(i),1),[2 3]), ...
        '13');
        H0=H0+Q+Q';
     end
  end

