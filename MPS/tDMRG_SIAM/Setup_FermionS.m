%testfile: setup generic model with spinfull fermions and NN-couplings and
%calculate ground state using standard DMRG

%BB 11-2016; QSpace3.0

global para

%basic model parameters
para.N=10;			     %length
para.t=1*ones(para.N,1);             %hopping
para.U=0*ones(para.N,1);             %on-site coulomb interaction strengh
para.V=0*ones(para.N,1);             %potential
para.hz=0*ones(para.N,1);            %magnetic field
para.mu=0*ones(para.N,1);            %chemical potential
para.model='FermionS_NN';               %model specification: spinful fermions
para.sym='Acharge,SU2spin';          %QSpace symmetries
para.flavor=1;			     %number of channels
para.Lambda=1;

%Symmetry checks
if (isequal(para.sym,'Acharge,SU2spin') | isequal(para.sym,'Pcharge,SU2spin') | isequal(para.sym,'SU2charge,SU2spin') | isequal(para.sym,'SU2spin')) && any(para.hz)
   error('ERR: SU2spin symmetry not compatible with finite hz');
elseif (isequal(para.sym,'SU2charge') | isequal(para.sym,'SU2charge,Aspin') | isequal(para.sym,'SU2charge,SU2spin')) && (any(para.V) || any(para.mu))
   error('ERR: SU2charge symmetry not compatible with finite chemical potential and/or onsite potential');
end

%DMRG parameters
para.update='2bond';                 %specifies DMRG update procedure - always '2bond' in case of symmetries!
para.NL=10;                          %dimension of Lanczos subspace
para.NRec=4;                         %maximum number of restarts for given optimization step
para.Ltol=1e-12;                     %tolerance of Lanczos alg
para.swmax=4;                        %maximum number of sweeps
para.sw=1;                           %sweep counter
para.precision=1e-8;                 %DMRG convergence tolerance
para.stol=1e-5;                      %minimum singular value kept in orthogonalization step
para.Nkeep=256;                      %max. bond dimension for non-symmetric calculation
para.eig_sigma='sa';                 %Calculate smallest ('sa') eigenvalue in optimization problem

%saving options
para.savedexist=0;
para.sflag=0;

disp(para);

%---------------------------------------------------------------------%
% CONVENTION: { L, R, sigma } order
% CONVENTION: first index of operators, ... connects to conj()
%---------------------------------------------------------------------%

%operator setup
[op.FC,op.Z,op.S,I]=getLocalSpace('FermionS',para.sym, 'NC',para.flavor);
for k=1:length(op.FC)
   op.N(k) = QSpace(contractQS(op.FC(k),'13*',op.FC(k),'13'));
   op.FC(k).info.itags{3} = 'm*';
end
op.ID=I.E;
op.Vac = getvac(op.Z,2);

%generate local Hamiltonian terms
op = GenHamiltonian(op);

%Operator storage for DMRG procedure
%Stores hamiltonian of left-right part
op.HLR=QSpace; op.HLR=op.HLR(ones(para.N+1,1));

%Stores next neighbour couplings
op.CLR=cell(para.N,length(op.h2l(1,:)));
for k=1:para.N+1
  for l=1:length(op.h2l(1,:))
    op.CLR{k,l}=QSpace;
  end
end


%perform simple NRG run to initialize MPS
H0 = op.h1{1};
H0.info.itags = {'','*'};
A0=QSpace(permuteQS(getIdentityQS(op.Vac,1,op.ID,1),[1 3 2]));

FF=op.FC;
NN=op.N;
ff=para.t(2:end);
if (isequal(para.sym,'Acharge,Aspin') | isequal(para.sym,'Pcharge,Aspin') | isequal(para.sym,'Aspin'))
    gg=zeros(para.N-1,3);
    gg(:,1)=para.mu(2:end)+para.V(2:end)-para.hz(2:end)./2;
    gg(:,2)=para.mu(2:end)+para.V(2:end)+para.hz(2:end)./2;
    gg(:,3)=para.U(2:end);
    NN(3) = op.N(1)*op.N(2);
    for i=1:numel(FF),
        if isequal(size(FF(i).Q,2),3), FF(i).Q(3)=[]; FF(i).info={}; end
    end
elseif (isequal(para.sym,'Acharge,SU2spin') | isequal(para.sym,'SU2spin') | isequal(para.sym,'Pcharge,SU2spin'))
        gg=zeros(para.N-1,2);
    gg(:,1)=para.mu(2:end)+para.V(2:end);
    gg(:,2)=para.U(2:end);
    NN(2) = (op.N(1)-op.ID)*op.N(1);
elseif (isequal(para.sym,'SU2charge,Aspin') | isequal(para.sym,'SU2charge'))
        gg=zeros(para.N-1,3);
    gg(:,1)=-para.hz(2:end)./2;
    gg(:,2)=+para.hz(2:end)./2;
    gg(:,3)=para.U(2:end);
    NN(3) = op.N(1)*op.N(2);
elseif isequal(para.sym,'SU2charge,SU2spin') 
        gg=zeros(para.N-1,1);

    NN(1) = (op.N(1)-op.ID)*op.N(1);
    gg(:,1)=para.U(2:end);
end    

Nkeep_NRG=para.Nkeep(1)*ones(para.N,1);  
opts_nrg={'Nkeep',para.Nkeep,'zflag',1};
[NRG,Inrg]=NRGWilsonQS(H0,A0,para.Lambda,ff,FF,op.Z,gg,NN,opts_nrg{:});

%Start with MPS ground state obtained from NRG run
[A,IA] = GetGS_NRG(NRG,Inrg);

clearvars -except A op para



for k=para.N-1:-1:1
    opts = {'stol',para.stol,'Nkeep',para.Nkeep,'nb',k};
    [A(k),A(k+1),SV]=MPSOrthoQS(A(k),A(k+1),'<<',opts{:});
    op=OpUpdate(A(k+1),op,k+1,'<<');
end


results = struct;
para.sw=1;
[A,op,results]=GS_DMRG(A,op,results,'-q');

op = TrotterGates(op,1i*0.025);
op.TGate2 = op.TGate; %for half steps in second order decomposition
op = TrotterGates(op,1i*0.05); %for full steps in second order decomposition


[A1,~,I] = TrotterEvol(A,op);
