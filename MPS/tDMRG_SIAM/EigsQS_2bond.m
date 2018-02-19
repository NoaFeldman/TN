function [B1,B2,Ie,SVD_Info]=EigsQS_2bond(B1,B2,op,k,lrdir,varargin)
%Two-site DMRG optimization using a standard Lanczos alg. in QSpace
%bond-update procedure from A.Weichselbaum Non-abelian symmetries in tensor networks:
%see "Prof. Dr. Peter Arbenz lecture note (ETHZ)" for detailed description of Lanczos implementation
%http://people.inf.ethz.ch/arbenz/ewp/lnotes.html
%Input:
%- B1,B2: describe neighbouring A matrices of the MPS
%- lrdir: sets direction for canonical update ('>>' for right-move, '<<' for left-move)
%- op: structure containing local Hamiltonian terms as well as global left and right Hamiltonian parts
%- k: current bond to be udpated
%- varargin: options to be passed
%            'Nkeep': number of kept state in SVD
%            'stol': minimum singular value kept in SVD
%            'NL': size of Krylov space
%            'NRec': nunmber of Lanczos restarts
%            'Ltol': convergence tolerance for Lanczos algorithm
%            '-q': quiet-mode - information of Lanczos optimization not printed
%Output:
%- B1,B2: optimized neighbouring A matrices of MPS in canonical form
%- Ie: ground-state energy calculated from local update
%- SVD_Info: Info-structure containing information about truncation (SVD spectrum ect.)

% BB 11-2012, 11-2016; QSpace 3.0

global  para

%flag: smallest (sa) or largest (la) ev
if isfield(para,'eig_sigma') para.eig_sigma='sa'; end
%number of Lanczos restarts
if isfield(para,'NRec') para.NRec=5; end
%eigenvalue tolerance
if isfield(para,'Ltol') para.Ltol=1e-10; end

%input options (see above)
getopt('INIT',varargin);
qflag  = getopt('-q'); 
if ~isempty(find(strcmp('stol', varargin)))  stol=getopt('stol',varargin); else stol=para.stol; end
if ~isempty(find(strcmp('Nkeep', varargin))) Nkeep=getopt('Nkeep',varargin); else Nkeep=para.Nkeep(1); end
if ~isempty(find(strcmp('NL', varargin)))    NL=getopt('NL',varargin); else Nkeep=para.NL; end
if ~isempty(find(strcmp('NRec', varargin)))  NRec=getopt('NRec',varargin); else NRec=para.NRec; end
if ~isempty(find(strcmp('Ltol', varargin)))  Ltol=getopt('Ltol',varargin); else Ltol=para.Ltol; end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Projection in bond picture (combining index Ls and sR)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

bond = 'A-A';
idx_c=[1,2]; %contraction index for orthogonalization

%bond: two site combination
A=QSpace;
A=QSpace(contractQS(B1,2,B2,1));%2site combination

%setting local indices
k1=k; k2=k+1;

%number of coupling terms
m=size(op.h2l,2);

%setting up identity matrices
%E1=QSpace(getIdentityQS(A,1,A,2));
E1=QSpace(getIdentityQS(A,1,op.h1{k1}));
E1.info.itags{3} = 'L*';
%E2=QSpace(getIdentityQS(A,3,A,4));
E2=QSpace(getIdentityQS(A,3,op.h1{k2}));
E2.info.itags{3} = 'R*';

%projecting 2-site block in effective basis
A=contractQS(conj(E1),A);
A=QSpace(contractQS(A,conj(E2)));


%projecting operators in effective basis
%left block
Q=QSpace; ic=[1,2];
if ~isequal(k1,1) && ~isempty(QSpace(op.HLR(k1-1)))
    op.HL_bond=QSpace(contractQS(conj(E1),ic, contractQS(op.HLR(k1-1),2,E1,1), ic));
else
    op.HL_bond=QSpace;
end

%left coupling to k1
if ~isequal(k1,1)
    for i=1:m
        Q=QSpace(contractQS(conj(E1),ic, ...
            contractQS(op.CLR{k1-1,i},contractQS(op.h2l{k1,i},E1)), ...
            ic));
        op.HL_bond=op.HL_bond+Q;
    end
end
%local term at k1
if ~isempty(op.h1{k1})
    Q=contractQS(E1,op.h1{k1},[1,3,2]);
    op.HL_bond=op.HL_bond+QSpace(contractQS(conj(E1),ic,Q,ic));
end

%right block
ic=[1,2]; ix=[1,3];
if ~isequal(k2,para.N) && ~isempty(QSpace(op.HLR(k2+1)))
    op.HR_bond=QSpace(contractQS(conj(E2),ic, contractQS(op.HLR(k2+1),E2), ic));
else
    op.HR_bond=QSpace;
end
%right coupling to k2
if ~isequal(k2,para.N)
    for i=1:m
        Q=QSpace(contractQS(conj(E2),ic, ...
            contractQS(op.CLR{k2+1,i}, contractQS(op.h2r{k2,i},E2)), ... % RLs
            ic));
        op.HR_bond=op.HR_bond+Q;
    end
end
%local term at k2
if ~isempty(op.h1{k2})
    Q=contractQS(E2,op.h1{k2},[1,3,2]);
    op.HR_bond=op.HR_bond+QSpace(contractQS(conj(E2),ic,Q,ic));
end

%intersite coupling
for i=1:m
    if ~isempty(op.h2r{k1,i})
        op.h2r_k1{i}=contractQS(conj(E1),ic, contractQS(E1,op.h2r{k1,i},[1,3,2,4]), ic);
    else
        op.h2r_k1{i}=QSpace;
    end
    if ~isempty(op.h2l{k2,i})
        op.h2l_k2{i}=contractQS(conj(E2),ic, contractQS(E2,op.h2l{k2,i},[1,3,2,4]), ic);
    else
        op.h2l_k2{i}=QSpace;
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Calculation of smallest (biggest) eigenvalue%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%starting vector
A1=A;
n=1;
while n <= para.NRec %repeat with new found ritz vector until convergence is reached
    
    clear alpha beta R T B W Q
    
    %Initialization of Lanczos steps
    Q={A1};
    %forming B = H*A
    if isequal(bond,'A-A')
        R=ContractHQS_2bond(A1,op);
    end
    
    %orthogonalization using Gram-Schmidt procedure
    a=contractQS(conj(A1),idx_c,R,idx_c);
    alpha(1)=a.data{:};
    if ~isscalar(a)
        fprintf('\n ERR <A|B> not a scalar: %.10g\t', a);
    end
    R=R-alpha(1)*A1;
    R=QSpace(R);
    b=contractQS(conj(R),idx_c,R,idx_c);
    if ~isempty(QSpace(b))
        beta(1)=sqrt(b.data{:});
    else
        Ie=alpha;
        break;
    end
    
    
    %check if MPS is already sufficiently
    if beta(1)<para.Ltol
        Ie=alpha;
        if ~qflag
            fprintf('\n abort Lanczos after step %g \t',1);
        end
        break;
    end
    
    if ~isscalar(b)
        fprintf('\n ERR  <B|B> not a scalar: %.10g\t', b);
    end
    
    
    
    %Iterative Building of Krylov space
    for j=2:para.NL
        W=Q{j-1}; B=QSpace(R)*(1/beta(j-1)); Q={Q{:},B};
        %build next Krylov space
        if isequal(bond,'A-A')
            R=ContractHQS_2bond(B,op);
        end
        R=R-beta(j-1)*W;
        a=contractQS(conj(B),idx_c,R,idx_c);
        alpha(j)=a.data{:};
        R=R-alpha(j)*B;
        
        %%%full orthogonalization (useful if Lanczos becomes unstable)
%         for l=1:j-1
%             a=contractQS(conj(Q{l}),idx_c,R,idx_c); if ~isempty(a.data)
%                 R=R-a.data{:}*Q{l};
%             end
%         end
        %%%
        
        R=QSpace(R);
        b=contractQS(conj(R),idx_c,R,idx_c);
        if ~isempty(QSpace(b))
            beta(j)=sqrt(b.data{:});
        else
            beta(j)=0;
        end
        %stop if weights get too small
        if beta(j)<para.Ltol
            if ~qflag
                fprintf('\n abort Lanczos after step %g of recursion %g; Cweight=%g <= Ltol=%g \n',j,n,beta(j),para.Ltol);
            end
            break;
        end
    end
    
    %calc of ritz value and ritz vector
    T=diag(beta(1:end-1),1);
    T=T+T'+diag(alpha);
    
    %options for matlab eigs
    opts.disp = 0;
    opts.tol = para.Ltol;
    opts.issym =1;
    sigma=para.eig_sigma;
    
    %solving for smallest eigenvalue
    [EV,Ie]=eigs(T,1,sigma,opts);
    
    %building ritz vector
    A1=QSpace;
    for i=1:size(Q,2)
        A1=A1+EV(i)*Q{i};
    end
    A1=QSpace(skipZerosQS(A1));
    
    %normalizing ritz vector
    norm=contractQS(A1,[1,2],conj(A1),[1,2]);
    A1=A1*(1./sqrt(norm.data{:}));
    
    %convergence check: decide whether to run the alg. again
    A2=QSpace;
    if isequal(bond,'A-A')
        A2=ContractHQS_2bond(A1,op);
    end
    
    res=A2-Ie*A1;
    delta=contractQS(conj(QSpace(res)),idx_c,QSpace(res),idx_c);
    delta=sqrt(delta.data{:});
    n=n+1;
    
    %if delta suffciently small stop else repeat algorithm with newly built A1
    if delta<para.Ltol 
        if ~qflag
           fprintf('\n abort Lanczos after recursion %g; Cweight=%g <= Ltol=%g \n',n,delta,para.Ltol);    
        end
        break;
    end
end

if n>para.NRec && ~qflag
    fprintf('\n Lanzcos did not converge in maximum number of recursions NRec=%g; Cweight=%g >= Ltol=%g \n',para.NRec,delta,para.Ltol);
end

%%%%%%%%%%%%%%%%%%%%%%%%
%truncation step       %
%%%%%%%%%%%%%%%%%%%%%%%%
opts = {'stol',para.stol,'Nkeep',para.Nkeep(1)};

%calculate SVD
[B2,B1,SVD_Info]=orthoQS(A1,1,lrdir,opts{:});
%relabel bond itag
if isequal(lrdir,'>>')
    B1.info.itags{1} =  ['B',num2str(k),'*'];
    B2.info.itags{1} =  ['B',num2str(k)];
elseif isequal(lrdir,'<<')
    B1.info.itags{1} =  ['B',num2str(k)];
    B2.info.itags{1} =  ['B',num2str(k),'*'];
end
B1=QSpace(contractQS(E1,B1, [1 3 2]));
B2=QSpace(contractQS(B2,E2));

end %function



function B=ContractHQS_2bond(A,op)
%Multiply Hamiltonian to local 2-site bond mps A

m=size(op.h2r,2);
Q=QSpace; B=QSpace;
%contract A to left and right hamiltonian and add up
B=QSpace(contractQS(op.HL_bond,A));
Q=QSpace(contractQS(A,op.HR_bond));
B=B+Q;

%intersite coupling
for i=1:m
    Q=contractQS(op.h2r_k1{i},A);
    Q=contractQS(Q,op.h2l_k2{i});
    B=B+QSpace(Q);
end
B=skipZerosQS(B);
end
