function [results] = Corr_tDMRG(Ket,op,corr,results,varargin)
%Determines time-dependent correlator < f1_{idx1}(t) f2_{idx2} > of
%fermionic chain model with NN-interactions only using standard tDMRG
%Input:
%- Ket: initial MPS to be time-evolved
%- op: structure containing local operators
%- corr: structure containing information about sites (corr.idx1,corr.idx2)
%  and operators (corr.f1,corr.f2) for time-dependent correlator
%  <f1_{idx1}(t) f2_{idx2}>; IMPORTANT in the presence of additional (third)
%  multiplet indices in f1,f2: convention is chosen such that third operator 
%  index is ALWAYS ingoing (outgoing) for f1 (f2); otherwise code will produce error! 
%- results: structure for main results of simulation
%- varargin: options to be passed
%            'Nkeep': number of kept state in SVD
%            'stol': minimum singular value kept in SVD
%            '-q': quiet-mode
%            'tau': Trotter time-step
%            'Tm': time after which measurement takes place (default Tm=tau)
%            'T': maximum time to which the time evolution is carried out
%Output:
%- A: optimized MPS
%- Z: normalization weigth (imaginary time-evolution)
%- Info: structure containing information about truncation ect

%BB 10-2015, 11-2016; QSpace3.0

global para

%input parameters
getopt('INIT',varargin);
qflag  = getopt('-q');
if ~isempty(find(strcmp('stol', varargin))) stol=getopt('stol',varargin); else stol=para.stol; end
if ~isempty(find(strcmp('Nkeep', varargin))) Nkeep=getopt('Nkeep',varargin); else Nkeep=para.Nkeep(1); end
if ~isempty(find(strcmp('tau', varargin))) Nkeep=getopt('tau',varargin); else tau=para.tau; end
if ~isempty(find(strcmp('Tm', varargin))) Tm=getopt('Tm',varargin); else Tm=tau; end
if ~isempty(find(strcmp('T', varargin))) Nkeep=getopt('T',varargin); else T=para.T; end

Braflag=0; %=1 time evolution for Bra layer carried out explicitly; =0 assumes ground state calculation and multiplies bra with exp(-i nTm E0) for each measurement

%number of time steps necessary
N = floor(T/Tm);
results.tCorr = nan*ones(N+1,1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% initialize time evolution at t=0 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%initialize Trotter gates
op = TrotterGates(op,1i*tau/2);
op.TGate2 = op.TGate; %for half steps in second order decomposition
op = TrotterGates(op,1i*tau); %for full steps in second order decomposition

Bra = Ket;

%apply f2 including corresponding Zs to Ket
for k=para.N:-1:corr.idx2+1
    Z = op.Z;
    Z.info.itags = {['s',num2str(k)],['s',num2str(k),'*']};
    Ket(k) = QSpace(contractQS(Ket(k),Z));
end
Ket(corr.idx2) = QSpace(contractQS(Ket(corr.idx2),corr.f2));

%restore right-canonical form of Ket
for k=1:1:para.N-1
    opts = {'stol',para.stol,'Nkeep',para.Nkeep,'nb',k};
    [Ket(k),Ket(k+1),I]=MPSOrthoQS(Ket(k),Ket(k+1),'>>',opts{:});
end
for k=para.N-1:-1:1
    opts = {'stol',para.stol,'Nkeep',para.Nkeep,'nb',k};
    [Ket(k),Ket(k+1),I]=MPSOrthoQS(Ket(k),Ket(k+1),'<<',opts{:});
end

%measure correlator at t=0
temp = Ket(corr.idx1:para.N);
for k=para.N:-1:corr.idx1+1
    Z = op.Z;
    Z.info.itags = {['s',num2str(k)],['s',num2str(k),'*']};
    Ket(k) = contractQS(Ket(k),Z);
end
Ket(corr.idx1) = contractQS(Ket(corr.idx1),corr.f1);
results.tCorr(1) = mpsOverlapQS(Ket,Bra);
Ket(corr.idx1:para.N) = temp;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% start tDMRG
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~qflag
    fprintf('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% \n ');
    fprintf('\n Starting time evolution \n');
    fprintf('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% \n ');
end

tic; m=1;
%truncation options
opts = {'stol',stol,'Nkeep',Nkeep,'Nt',round(Tm/tau)};

%carry out N time steps
for n=1:N
    
    if Braflag
        %time evolve bra one time-step
        [Bra,~,I] = TrotterEvol(Bra,op,opts{:});
        
        %storing MPS info
        results.Info_bra.svdtr(n)=sum(I.svdtr);
        results.Info_bra.Nkeep(n)=max(I.Nkeep);
        
        if ~qflag
            fprintf(' Maximum bond dimension: %g \n',max(I.Nkeep));
            fprintf(' Cumulated discarded weight: %g \n', sum(I.svdtr));
            fprintf('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% \n ');
        end
    end
    
    %time evolve ket by single time-step tau
    [Ket,~,I] = TrotterEvol(Ket,op,opts{:});
    
    %storing MPS info
    results.Info_ket.svdtr(n)=sum(I.svdtr);
    results.Info_ket.Nkeep(n)=max(I.Nkeep);
    
    if ~qflag
        fprintf(' Maximum bond dimension: %g \n',max(I.Nkeep));
        fprintf(' Cumulated discarded weight: %g \n', sum(I.svdtr));
        fprintf(' Time  t=%g \n', n*Tm);
        fprintf('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% \n ');
    end
    
    

    
    
    %measure correlator at t=n*tau
    temp_ket = Ket(corr.idx1:para.N);
    temp_bra = Bra(1);
    %apply phase factor
    if ~Braflag
        Bra(1) = Bra(1)*exp(-1i*n*Tm*para.E0);
    end
    for k=para.N:-1:corr.idx1+1
        Z = op.Z;
        Z.info.itags = {['s',num2str(k)],['s',num2str(k),'*']};
        Ket(k) = contractQS(Ket(k),Z);
    end
    Ket(corr.idx1) = contractQS(Ket(corr.idx1),corr.f1);
    results.tCorr(n+1) = mpsOverlapQS(Ket,Bra);
    Ket(corr.idx1:para.N) = temp_ket;    
    Bra(1) = temp_bra;
    results.t(n+1) = n*Tm;
    results.t_CPU(n+1) = toc;
    
    if isequal(para.sflag,1) && isequal(mod(n,30),0)
        save(para.filename,'para','results');
    end
    
    %abort if MPS reaches maximum number of states
    if n>5 && all(results.Info_ket.Nkeep(end-4:end)==Nkeep)
        save(para.filename,'para','results');
        fprintf('Abort real-time evolution because maximum bond dimension D=%g was reached after T=%g \n',Nkeep,results.t(n+1));
        break;
    end
    
end %time evolution for

end %function
