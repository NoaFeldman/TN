% -------------------------------------------------------------------- %
% setup operator set for standard Kondo problem (single channel)
% -------------------------------------------------------------------- %

% Note on coupling
% in SIAM, the coupling is given by sqrt(rho)*V =: sqrt(Gamma/pi)
% sqrt(rho) since only a single c_k is coupling to the dot
%
% Kondo coupling: symmetric SIAM neglecting excitation energies much
% smaller than U and U+epsd, one has J = 4V²/U resulting in an overall
% Kondo coupling of 8*Gamma/U between the impurity and the first
% site in the Wilson chain and regular exponentially decaying coupling
% thereafter (see handnotes)
%
% J_calc = (2*rho*J) * D

  global param

  if ~exist('J','var') || ~exist('B','var') || ~exist('Lambda','var')
     wblog('TST Setting test parameters.');
     B=0.01; J=0.03; setParam
  end

  param=struct( ...
     'fJ',J, 'B', B, 'Lambda', Lambda, 'N', N  ...
  );

% pg_r = pseudo-gag coefficient r
  co={};
  if exist('z',   'var'), co(end+1:end+2)={'z',z   }; param.z=z;    end
  if exist('pg_r','var'), co(end+1:end+2)={'r',pg_r}; param.r=pg_r; end

  [ff,q,Ic]=getNRGcoupling(1,Lambda,N,co{:});

  if exist('J','var')
       ff(1)=J;
  else ff(1)=8*Gamma/(pi*U); end

  param.fJ=ff(1);

  if exist('TK0')==1, TK=TK0; else
  TK=TKondo; end

  wblog('<i>', 'Kondo parameters (TK=%.4e)\N', TK);
  disp(param);

% -------------------------------------------------------------------- %
% operator setup
% -------------------------------------------------------------------- %

  initdef('QTYPE','QSz');
  initdef('NC',    1);   % number of channels
  initdef('ISP',NC/2);   % default Kondo impurity spin = 1/2

  m2=2*ISP+1;

  if exist('B','var')
     % like for pauli Sz, +mz ... -mz (pos. mz come first; see below)
       data = B * (+ISP:-1:-ISP);  % assume H_B=+B·Sz
     % data = B * (-ISP:ISP);
  else data = zeros(1,m2); end

% spin matrizes
  s0 = spinmat( 2*ISP+1); % regular Sz (incl. prefactors)
  sx = spinmat('-pauli' ); % plain Pauli matrizes

% make real assuming only sy*sy will appear (!)
  s0{2}=+1i*s0{2};
  sx{2}=-1i*sx{2};

% -------------------------------------------------------------------- %
  if isequal(QTYPE,'QSz')
% conservation of number of particles (charge Q) and Sz_tot
% NB! particle number in each channel is conserved
% -------------------------------------------------------------------- %

  % Hamiltonian H0 and its basis A0
  %    [ charge, 2*Sz ] basis // up down
    QS = [-1  0              %    0  0
           0 +1              %    1  0  % corresponds to standard
           0 -1              %    0  1  % Pauli Sz convention
          +1  0 ];           %    1  1

    Q0 = [ zeros(m2,NC), 2*(ISP:-1:-ISP)' ]; % ie. 1,-1 for spin 1/2

    Qloc=repmat({zeros(4,NC+1)},1,NC);
    for i=1:NC, Qloc{i}(:,[i end])=QS; end

  elseif isnumeric(QTYPE) && isscalar(QTYPE)
    error('Wb:KONDO','not implemented for QTYPE=%d',QTYPE);
  else
    QTYPE % display
    error('Wb:SIAM','Invalid QTYPE');
  end

  cs=QSpace(2,NC); cn=QSpace(2,NC); cz=QSpace(1,NC);
  for l=1:NC
     Ql=Qloc{l};
     cs(:,l)=[ QSpace(...
             Ql([1 2],:),  1, ...  % destroy spin-up: d(Sz)=-1
             Ql([3 4],:),  1  ...
          )
          QSpace(...
             Ql([1 3],:),  1, ...  % destroy spin-down: d(Sz)=+1
             Ql([2 4],:), -1  ...
          )
     ];

     cn(:,l)=[ QSpace(Ql([3 4],:),  1)      % c_up * n_down
          QSpace(Ql([2 4],:), -1) ];   % c_down * n_up

     zs(l)=QSpace(Ql, diag((-1).^[0 1 1 2]), 'operator');
  end

% -------------------------------------------------------------------- %
% -------------------------------------------------------------------- %

  H0=QSpace( ...
     { Q0, Q0 }, ...
     mat2cell( data, 1, ones(m2,1) )' ...
  );

  [A0,Aloc]=QSpace(Q0,Qloc{:},'identity');
  E0=QSpace(Q0, eye(m2), 'operator'); % identity for Kondo spin

% E4=QSpace(QS, eye(4), 'operator');
  E4=contractQS(Aloc,1:NC,Aloc,1:NC);

  S0=repmat(QSpace,3,1); SX=repmat(QSpace,3,NC);
  for p=1:3
   % s0 => S0 include (1/2) factors
     S0(p)=QSpace(Q0,s0{p},'operator');
   % sx => SX DO NOT (1/2) factors (based on Pauli matrices)
   % => H_Kondo = 2*ff(1) * S.S
     for l=1:NC
     SX(p,l)=QSpace(Qloc{l}(2:3,:),sx{p},'operator'); end
  end

  ic=1:(2+NC); ic(2)=[]; % ic=[1 3];
  H0=QSpace(contractQS(A0,ic, contractQS(H0,2,A0,1),ic,'conjA'));

% see definition of TK in terms of J in TKondo.m
  Q=QSpace; f=ff(1);
  for l=1:NC, k=2+l;
     ix=cmat_iperm(2+NC,2+l); ix(2)=[];
     for p=1:3
        Q=Q+contractQS(A0, ic, ...
        contractQS(contractQS(f*S0(p),2,A0,1), 2+l, SX(p,l),2), ix,'conjA');
     end
  end
  H0=H0+Q;

% compactify local state space
  A0=QSpace(contractQS(A0,3:2+NC,Aloc,1:NC));

  ff=ff(2:end); % done with the first iteration!!

% -------------------------------------------------------------------- %
% general impurity / Wilson site operators
% -------------------------------------------------------------------- %

% annihilation operators in LOCAL state space (set of Wilson site)
  FC=repmat(QSpace,2,NC); Q=Aloc;
  for k=1:NC,
     iperm=cmat_iperm(NC+1,k); % rank(Aloc) = NC+1 with R to the right
     for i=1:2 % spin up/down
         FC(i,k)=...
           contractQS(Aloc, 1:NC,...
           contractQS(Q,k,cs(i,k),2), iperm(1:end-1),'conjA'...
         );
     end
     Q=contractQS(Q,k,zs(k),2,iperm);
  end

% Z operator for Wilson site (NRG)
  Z=contractQS(Aloc, 1:NC, Q, 1:NC,'conjA');

% --------------------------------------------------------------------
% operators for dmNRG
% -------------------------------------------------------------------- %

% Z operator for impurity ( dmNRG -> Z0*op[12] )
% tensorprod Z0 required for tensorprod op[12] below
  Z0=mpsTensorProdQS(E0,Z);
  wblog('CHK','\N(Z,E4) or (E4,Z) ??\N');

% NB! SX([12]) won't preserve symmetries -> dmNRG won't like it :/
% F=[]; C=SX(3); op2=SX(3);

% NB! extra factor of 2, e.g. in case of spin 1/2:
% S_ = SX-i*SY = QSpace([ 0 -1; 0 +1], 2);

  S_=S0(1)-S0(2); % see def. of S0(s0) above: s0{2}=+1i*s0{2};
  cc = [
     S_   ,  FC(2) % only for first channel
     S0(3), +FC(1)
     S_'  ,  FC(1)
     S0(3), -FC(2)
  ];

  m=size(cc,1); o=repmat(QSpace,m,1);
  for i=1:m
  o(i)=mpsTensorProdQS(cc(i,1),cc(i,2)); end

  o(5)=o(1)+o(2); o(6)=o(3)+o(4);

  op1=[];
  op2=param.fJ*o([ 5 6 ]); zflags=ones(size(op2));

  OP1=op1;
  for i=1:length(OP1)
     OP1(i)=contractQS(A0,[1 3],...
     contractQS(OP1(i),[3 4],A0,[1 3]),[1 2],'conjA');
  end

  OP2=op2;
  for i=1:length(OP2)
     OP2(i)=contractQS(A0,[1 3],...
     contractQS(OP2(i),[3 4],A0,[1 3]),[1 2],'conjA');
  end

  FC_=FC; % rdma load FC so it will be overwritten

% -------------------------------------------------------------------- %
  clear k i p ic ix f m2 Q % QS
% -------------------------------------------------------------------- %

