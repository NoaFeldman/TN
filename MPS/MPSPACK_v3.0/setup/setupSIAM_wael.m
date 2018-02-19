% -------------------------------------------------------------------- %
% setup general parameters
%  Mit diesem Programm wird der Hamiltonian H=H_free+H_tunneling+H_Hund beschrieben.
%  Wobei H_Hund=-JH*S_d*S_d.
%  Die Idee ist die Struktur von setupKondo zu uebernehmen und statt mehrere
%  Leitungsbaender (kanaele NC) zu betrachten, betrachten wir mehrere Impurity-Levels (NL).
%  Es wird also das Huepfen von einem Leitungsband auf NL-Impurity levels (Dots) beschrieben.
% -------------------------------------------------------------------- %

  global param

  if ~exist('U','var') && ~exist('Lambda','var') % U=0.12; epsd=-U/2;
     Gamma=0.01; Lambda=2; N=50; B=0; Nkeep=512;
     JH=0.000075; NL=3; Dz=0;
  end
  initdef('B',0);

% TEST
% JH=0; NL=1; Gamma=0.01; Nkeep=512; N=50;
  JH=0.02; NL=3; Gamma=0.01; Nkeep=512; N=90; Lambda=2;

  param=struct( ... % 'U', U, 'epsd', epsd, 
      'Gamma', Gamma, 'B', B, 'JH', JH , 'NL', NL, 'Dz', Dz, ...
      'Nkeep', Nkeep,'N', N, 'Lambda', Lambda  ...
  );

  if isstr(B), B=eval(B); param.B=B; end

  co={};
  if exist('z',   'var'), co(end+1:end+2)={'z',z   }; param.z=z;    end
  if exist('pg_r','var'), co(end+1:end+2)={'r',pg_r}; param.r=pg_r; end

  ff=getNRGcoupling(Gamma,Lambda,N,co{:});
  
   s0 = spinmat(2);
	s0{2}=-1i*s0{2};
  
  TK=exp(-JH/Gamma); % estimate
  global g_TKondo ; g_TKondo=TK;

  wblog('<i>', 'parameters for SIAM/Anderson hybrid (TK ~ %.4g)\N',TK);
  disp(param);

% -------------------------------------------------------------------- %
% operator setup
% -------------------------------------------------------------------- %

   % Hamiltonian H0 and its basis A0
   %    [ charge, 2*Sz ] basis // up down
     QS = [-1  0              %    0  0    empty
            0 -1              %    0  1    spin down
            0 +1              %    1  0    spin up
           +1  0 ];           %    1  1    doubly occupied

% -------------------------------------------------------------------- %
% NEW-Begin!!!
% -------------------------------------------------------------------- %


% site annihilation operators
  FC(1)=QSpace(...  % spin up
          [-1  0; 0  1],  1, ...
          [ 0 -1; 1  0],  1  ...
        );

  FC(2)=QSpace(...  % spin down
          [-1  0; 0 -1],  1, ...
          [ 0  1; 1  0], -1  ...
        );

 % Z-Operator!!!
   Z = QSpace(...
          { QS, QS }, ...
          mat2cell( (-1).^[0 1 1 2], 1, ones(4,1))' ...
       );


% -------------------------------------------------------------------- %
% -------------------------------------------------------------------- %

  if NL>1
     Q=repmat({QS},1,NL);
     A=QSpace(Q{:}, '-Rlast', 'identity');
  elseif NL==1
     A=QSpace(QS,eye(size(QS,1)),'operator');
  else
  error('Wb:ERR','invalid NL=%g',NL); end

% Die Spin-Operatoren fuer die drei Levels in der A-Basis :
% SX for the three levels :
  S0=QSpace(3,NL); ic=1:NL; clear pp

  for p=1:3; s0{p}=QSpace(QS(2:3,:),s0{p},'operator'); 
  for l=1:NL;
     pp([ l, 1:l-1, l+1:NL ])=1:NL;
     S0(p,l)=QSpace(contractQS(A,ic, ...
     contractQS(s0{p},2,A,l), pp, 'conjA'));
  end
  end

  s0=cat(1,s0{:}); SX=S0(:,1);

  for p=1:3
  for l=2:NL, SX(p,1)=SX(p,1)+S0(p,l); end, end

  SS=QSpace(SX(1))*SX(1)-SX(2)*SX(2)+SX(3)*SX(3);
 
% the Dot-Hamiltonian with JH the Hund coupling : (64x64)
  Hdot=-JH*SS + Dz*SX(3)*SX(3);
 
% dZ in der A-Basis mit Kontraktion von den Z-Operatoren, in der Konv. f0d3d2d1|0>!!!
  FD=QSpace(2,NL);
  for l=1:NL, pp([NL:-1:l, 1:(l-1)])=1:NL;
  for i=1:2 % spin
     Q=contractQS(FC(i),2,A,l);
     for j=l+1:NL, Q=contractQS(Z,2,Q,j); end
     FD(i,l)=QSpace(contractQS(A,ic,Q,pp,'conjA'));
  end
  end

  % if U~=0
  %    ND=QSpace; % Gesamtteilchenzahl am Dot
  %    for i=1:numel(FD), ND=ND+FD(i)'*FD(i); end
  %    Hdot=Hdot+(U/2)*ND*ND;
  % end

% hier kommt die Aufbauschleife fuer Q!
  A0 = QSpace(A.Q{end}, QS, 'identity');

% Hdot in der A0-Basis 256x256 :
  H0=contractQS(A0, [1 3], contractQS(Hdot, 2, A0, 1), [1 3], 'conjA');
  H0=0*QSpace(contractQS(A0,[1 3],A0,[1 3])) + H0; % keep this!
 
% Addition des Dot-Hamiltonians!!! 
% Definition des Huepf-Hamiltonians :
% d*f0 in der A0 Basis fuer die verschiedenen Dot-Levels!!!! 
  for i=1:2;
  for l=1:NL;
  % Multiplikation mit der Kopplung :
    Q=ff(1)*QSpace(contractQS(A0, [1 3], ...
    contractQS(contractQS(FD(i,l),1,A0,1),3,Z*FC(i),2),[1 3], 'conjA'));
    H0=H0+Q+Q'; % i.e. d^{dag}*c+h.c.
  end
  end

% erster Huepfterm ist erledigt!
  f1=ff(1); ff=ff(2:end);

  NKEEP=repmat(2*Nkeep,1,5);

% Z Operator for dmNRG !!!!!!!

  Q=A; for i=NL:-1:1, Q=contractQS(Z,2,Q,NL); end % XXX A->Q XXX
  Z0=QSpace(contractQS(A,ic,Q,ic,'conjA'));

  op1=[]; op2=FD(:); % nostore=1; % => rdma.m

  wnrg=''; clear QS Q i p l

% -------------------------------------------------------------------- %
% -------------------------------------------------------------------- %

