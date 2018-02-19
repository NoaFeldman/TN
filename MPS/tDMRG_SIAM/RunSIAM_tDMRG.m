function RunSIAM_tDMRG(N,U,epsd,Gamma,hz)
%calculates spectral function of SIAM using tDMRG
%Input:
%- N: length of Wilson chain (counts all sites including impurity!)
%- U: coloumb interaction strength on impurity
%- epsd: onsite-energy on impurity
%- Gamma: hybridisation strength
%- hz: magnetic field on impurity

%BB 11-2016; QSpace 3.0

if isdeployed           % takes care of command line arguments
    if ischar(U), U = str2num(U); end
    if ischar(epsd), epsd = str2num(epsd); end
    if ischar(Gamma), Gamma = str2num(Gamma); end
    if ischar(hz), hz = str2num(hz); end
    if ischar(N), N = str2num(N); end
end

global para

%get Wilson chain parameters (example)
para.N=N;
para.Lambda=1.2;
para.Gamma = Gamma;
co={'-x','-AL'};
[ff,N,Ic]=getNRGcoupling(para.Gamma,para.Lambda,para.N+1,co{:});
para.t=ff;

%basic model parameters
para.U0=U;
para.U=0*ones(para.N,1);             %on-site coulomb interaction strengh (only impurity site!)
para.U(1)=para.U0;
para.V=0*ones(para.N,1);             %potential
para.hz0=hz;
para.hz=0*ones(para.N,1);            %magnetic field (only impurity site!)
para.hz(1)=para.hz0;
para.mu=0*ones(para.N,1);            %chemical potential (only impurity site!)
para.epsd = epsd;
para.mu(1)=-epsd;
para.model='FermionS_NN';            %model specification: spinful fermions
para.sym='Acharge,SU2spin';          %QSpace symmetries
para.flavor=1;                       %number of channels

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
para.swmax=10;                       %maximum number of sweeps
para.sw=1;                           %sweep counter
para.precision=1e-9;                 %DMRG convergence tolerance
para.stol=1e-5;                      %minimum singular value kept in orthogonalization step
para.Nkeep=512;                      %max. bond dimension for non-symmetric calculation
para.eig_sigma='sa';                 %Calculate smallest ('sa') eigenvalue in optimization problem

%tDMRG parameters
para.Nkeep_tDMRG=1024;
para.stol_tDMRG=1e-4;
para.tau=0.05;                       %Trotter time step
para.Tm=0.25;                         %measurement times
para.T=80;                            %maximum time for tDMRG evolution


%saving options & folder
para.sflag=1;
para.folder=strcat('results_SIAM/tDMRG/',sprintf('/N%dU%.5gGamma%.5ghz%.5gT%.5gtau%.5g',para.N,para.U(1),para.Gamma,para.hz(1),para.T,para.tau));
para.filename=strcat(para.folder,'/results.mat');

if ~exist(para.filename,'file')
    mkdir(para.folder);
end

user_name =getenv('USER');
para.RC_STORE=strcat(['/data/',user_name,'/'],para.folder);
para.sflag=1;

% //Set Enviroment
setenv('RC_STORE', para.RC_STORE);
setenv('CG_VERBOSE', '0');
if ~exist(para.RC_STORE,'dir')
    mkdir(para.RC_STORE);
else
    rmdir(para.RC_STORE,'s');
    mkdir(para.RC_STORE);
end

disp(para);


%---------------------------------------------------------------------%
% CONVENTION: { L, R, sigma } order
% CONVENTION: first index of operators, ... connects to conj()
%---------------------------------------------------------------------%

%operator setupS(
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



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GROUND STATE DMRG %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%bring MPS in right-canonical form
for k=para.N-1:-1:1
    opts = {'stol',para.stol,'Nkeep',para.Nkeep,'nb',k};
    [A(k),A(k+1),SV]=MPSOrthoQS(A(k),A(k+1),'<<',opts{:});
    op=OpUpdate(A(k+1),op,k+1,'<<');
end


results = struct;
para.sw=1;

%calculate ground state using DMRG
[A,op,results]=GS_DMRG(A,op,results,'-q');
para.E0 = results.E(end);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% tDMRG %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%setup correlator structure to calculate < f_1(t) f'_1 > using Corr_tDMRG.m
%Note: "corr" structure containing information about sites (corr.idx1,corr.idx2)
%  and operators (corr.f1,corr.f2) for time-dependent correlator
%  <f1_{idx1}(t) f2_{idx2}>; IMPORTANT in the presence of additional (third)
%  multiplet indices in f1,f2: convention is chosen such that third operator
%  index is ALWAYS ingoing (outgoing) for f1 (f2); otherwise code will produce error!

para.Nkeep=para.Nkeep_tDMRG;
para.stol =para.stol_tDMRG;

%setup correlator structure to calculate < f_1(t) f'_1 >
results_ffdag=struct;
corr.f1 = op.FC(1);  corr.f1.info.itags = {'s1','s1*','m*'};
corr.f2 = op.FC(1)';   corr.f2.info.itags = {'s1','s1*','m'};
corr.idx1=1; corr.idx2=1;
[results_ffdag] = Corr_tDMRG(A,op,corr,results_ffdag,'Tm',para.Tm);

%setup correlator structure to calculate < f'_1(t) f_1 >
results_fdagf=struct;
corr.f1 = op.FC(1)';  corr.f1.info.itags = {'s1','s1*','m'};
corr.f2 = op.FC(1);   corr.f2.info.itags = {'s1','s1*','m*'};
corr.idx1=1; corr.idx2=1;
J1 = QSpace(getIdentityQS(corr.f1,3,'-0'));
corr.f1 = QSpace(contractQS(corr.f1,3,J1,'1*'));
corr.f2 = QSpace(contractQS(corr.f2,3,J1,'1'));
[results_fdagf] = Corr_tDMRG(A,op,corr,results_fdagf,'Tm',para.Tm);



%calculate spectral function
GFp = results_fdagf.tCorr;
GFm = results_ffdag.tCorr;


if length(GFp)<length(GFm)
    GFm=GFm(1:length(GFp));
elseif length(GFp)>length(GFm)
    GFp=GFp(1:length(GFm));
end
GF = -1i*(GFp+conj(GFm));

%linear prediction
para.Nmax=10000;
[GF_LP] = LinearPrediction(GF,para.Nmax); GF_LP = GF_LP.';

t_LP=0:r.para.tau:Nmax*r.para.tau-r.para.tau;%
A_LP = trapz(t_LP,ones(length(omega),1)*GF_LP.*exp(1i.*omega'*t_LP),2);
A_LP = -(1/pi)*imag(A_LP);


results.GF = GF;
results.GFp=GFp;
results.GFm=GFm;
results.GF_LP = GF_LP;
results.A_LP = A_LP;

save(para.filename,'para','results');


end