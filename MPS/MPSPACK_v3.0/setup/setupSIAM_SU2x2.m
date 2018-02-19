% -------------------------------------------------------------------- %
% setup operator set for standard SIAM problem (single channel)
% using particle-hole and spin SU(2) symmetry (if B=0)
% Wb,Sep19,11
% -------------------------------------------------------------------- %

% adapted from setupSIAM.m

  global param

  setdef('U',0.08, 'Gamma',0.01, 'Lambda',2, 'N',50);
  setdef('epsd',-U/2, 'B',0);
  Bflag=(B~=0);

  istr=''; sym='';
  param=add2struct('-',istr,sym,U,epsd,Gamma,B,Lambda,N);

% pg_r = pseudo-gag coefficient r
  co={'-x'};
  if exist('z','var'), co(end+1:end+2)={'z',z}; param.z=z; end
  if exist('ALambda','var') && ~ALambda
     co(end+1)={'-w'}; % param.ALambda=0;
     if isfield(param,'ALambda') param=rmfield(param,'ALambda'); end
  end

  if exist('TK0')==1, TK=TK0; else
  TK=TKondo; end

  param.istr='one-channel Anderson model (SIAM)';

% -------------------------------------------------------------------- %
% operator setup
% -------------------------------------------------------------------- %

  o={'NC',1,'-v'}; clear ZFLAG*
  if epsd~=-U/2 && Bflag || isset('SYM') && isequal(SYM,'A,A')
     [FF,Z,SS,IS]=getLocalSpace('FermionS','Acharge,Aspin', o{:});

     NF=FF; % F3=FF;
     for i=1:numel(FF)
      % FF(i).Q=FF(i).Q(1:2); FF(i).info={};
        FF(i)=fixScalarOp(FF(i));
        NF(i)=FF(i)'*FF(i);
     end

     HU=skipzeros(NF(1)*NF(2)); % = Nup*Ndown
     H0=skipzeros(epsd*(NF(1)+NF(2)) + U*HU - B*SS(3)); % NB! want <Sz> >= 0!
     E0=IS.E; % QSpace(getIdentityQS(Z,1));

  elseif epsd~=-U/2 && ~Bflag || isset('SYM') && isequal(SYM,'A,SU2')
     [FF,Z,SS,IS]=getLocalSpace('FermionS','Acharge,SU2spin', o{:});
     E0=IS.E; % QSpace(getIdentityQS(Z,1));

     NF=QSpace(contractQS(FF,[1 3],FF,[1 3],'conjA'));
     HU=0.5*skipzeros(NF*NF-NF); % = (2/2)*Nup*Ndown
     H0=skipzeros(epsd*NF + U*HU);

  elseif epsd==-U/2 && Bflag || isset('SYM') && isequal(SYM,'SU2,A')
     ZFLAG=2;
     [FF,Z,SS,IS]=getLocalSpace('FermionS','SU2charge,Aspin', o{:});
     E0=IS.E; % QSpace(getIdentityQS(Z,1));

   % NB! 0.5*Q2 = (1/2)*(Nup+Ndown-1)^2 = (2/2)*(Nup*Ndown) - (1/2)(Nup+Ndown) + (1/2)
   % where the 2nd term already includes the onsite energy(!)
   % and the last terms is corrected by subtracting -E0/2
     HU=0.5*skipzeros(IS.Q2-E0);
     H0=skipzeros(U*HU-B*SS(3)); % NB! want <Sz> >= 0!

  else, ZFLAG=2;
     [FF,Z,SS,IS]=getLocalSpace('FermionS','SU2charge,SU2spin',o{:});
     E0=IS.E; % QSpace(getIdentityQS(Z,1));

   % NB! E0 simply corrects constant (irrelevant) energy offset
     HU=0.5*skipzeros(IS.Q2-E0); % = (2/2)*Nup*Ndown - (1/2) NF with the last
     H0=skipzeros(U*HU);         % already including the onsite energy
  end

  if isempty(H0), H0=0*E0; end

  FN=comm(FF,HU);

  if ~gotCGS(FF(1))
     A0=QSpace(zeros(1,size(Z.Q{1},2)), Z.Q{1}, 'identity');
  else
   % n=numel(FF); q=getsub(Z,1); for i=1:2, q.Q{i}(:)=0; end
     q=getvac(Z);
     A0=QSpace(permuteQS(getIdentityQS(q,1,Z,1),[1 3 2]));
  end

% A1=QSpace(getIdentityQS(A0,2,IS.E)); A1=permute(A1,[1 3 2]);
% A2=QSpace(getIdentityQS(A1,2,IS.E)); A2=permute(A2,[1 3 2]);

% ensure that subsequent calls use the same symmetries
% e.g. when run in a loop -- Wb,Oct21,11
  SYM=getsym(H0); param.sym=IS.sym;

  [ff,q,Ic]=getNRGcoupling(Gamma,Lambda,N,co{:});

% -------------------------------------------------------------------- %
  FC=FF; op1=[];

  if Bflag
     op1=[ FC(1), FN(1) ];
     op2=[ FC([1,1]) ]; zflags=ones(1,numel(FC));
  else
     op1=[ FC(1), FN(1), SS(end) ]; % SS(end): see SYM above!
     op2=[ FC(1), FC(1), SS(end) ]; zflags=[1 1 0];
  end

  if isset('ZFLAG') && ZFLAG>1
   % NB! FN uses HU which in case of particle-hole symmetry is 0.5*(N-E)^2!
   %    1/(omega - epsd - Delta - U*F/G)          % need to cancel epsd
   % => 1/(omega + U/2 - Delta - U*(F+0.5*G)/G)   % use F -> F+G/2
     ac_cmd='[ac,Ic]=getGC_bulla(ox, ax(:,1),  ax(:,2) + 0.5*(ax(:,1)),''-iB'');';
  else
     ac_cmd='[ac,Ic]=getGC_bulla(ox, ax(:,1),  ax(:,2),''-iB'');';
  end

  Z0=Z; % required for rdma.m // Wb,Sep21,13

% -------------------------------------------------------------------- %
% add first n0add Wilson sites to H0 and A0 // Wb,Jul01,13
% -------------------------------------------------------------------- %

  if isset('n0add') && n0add>0
     
     ZF=FF;
     for i=1:numel(FF), ZF(i)=contractQS(Z,2,FF(i),1); end

     for it=1:n0add
        wblog(' * ','adding %g/%g Wilson site to A0 ...',it,n0add);

        if isset('ZFLAG') && ZFLAG>1
         % NB! last site added must be opposite to ZFLAG used in NRGWilsonQS
           if mod(it+n0add-1,2)==(ZFLAG-2), q1=FF; else q1=ZF; end
           q2=q1;
        else
           q1=FF; q2=ZF;
        end

        for i=1:numel(q1)
           q1(i)=contractQS(A0,[1 3],contractQS(A0,3,q1(i),2),[1 3]);
        end

        A0_=A0; % remember last A0
        A0=QSpace(permuteQS(getIdentityQS(A0,2,IS.E,1),[1 3 2]));

      % propagate H0
        H0_=H0; % remember last Hamiltonian asscoiated with L of A0 now 
        H0=QSpace(contractQS(A0,[1 3],contractQS(H0,2,A0,1),[1 3]));

        for i=1:numel(q2)
           q1(i)=contractQS(q1(i),1,A0,1,'-L');  % F :: L(o)Rs
           q2(i)=contractQS(A0,3, ff(it)*q2(i),1,'-L'); % f' :: LRs(o)
           if numel(q2(i).Q)==2
                q=QSpace(contractQS(q1(i),[1 3],q2(i),[1 3]));
           else q=QSpace(contractQS(q1(i),[1 2 4],q2(i),[1 4 3]));
           end

           if isset('ZFLAG') && ZFLAG>1
                H0=H0+q;
           else H0=H0+q+q'; end
        end
     end

     fprintf(1,'\n'); clear ZF q1 q2
     param.n0add=n0add;

     ff=ff(n0add+1:end);
  end

% -------------------------------------------------------------------- %
% -------------------------------------------------------------------- %

  inl(1); disp(param);

