function [Phi,IPhi] = getGS_NRG(nrg_basis,nrg_info)
% usage:
% [Phi] = get_ground_state(nrg_info,nrg_basis)
%
% INPUT:
% nrg_info - info from NRG run
% nrg_basis - NRG basis
%
% OUTPUT:
% Phi - the ground state of the system
%
% Jun 2011, IW
% March 2015, BB - adapted to QSpaceV3

global para


N = nrg_info.N;

block_idx = 1; state_idx = 1;
[x,state_idx] = min(real(nrg_basis(N).HT.data{1}));
engs = nrg_basis(N).HT.data{1}(state_idx);
for n=2:length(nrg_basis(N).HT.data)
    if min(real(nrg_basis(N).HT.data{n})) < x,
        [x,state_idx] = min(real(nrg_basis(N).HT.data{n}));
        block_idx = n;
        engs = nrg_basis(N).HT.data{n}(state_idx);
    end
end
vect(1:length(nrg_basis(N).HT.data{block_idx})) = 0;
vect(state_idx) = 1;
% check if the energies were also sorted wrt to real part in NRG run
if state_idx ~= 1,
    printf('Warning: The energies during NRG runwere not sorted wrt real part');
end
% the last block of MPS related to the ground state
%AG = QSpace( [nrg_basis(N).HT.Q{1}(block_idx,:); nrg_basis(N).HT.Q{1}(block_idx,:)], vect.');
Q=QSpace( nrg_basis(N).HT);
Q=Q*0;
Q.data{block_idx}(state_idx)=1;

%Setup GS MPS
phi = QSpace(1,N);
for n=1:N-1
    phi(n) = QSpace(nrg_basis(n).AK);
    phi(n).info.itags = { strcat('B',num2str(n-1)),  strcat('B',num2str(n),'*') , strcat('s',num2str(n))};
end
phi(N) = QSpace(contractQS(nrg_basis(N).AT,2,Q,'2;*')); %
phi(N)=permuteQS(phi(N), [1 3 2]); % Indizes so permutieren, dass wieder LRs


%add vacuum state to MPS and turn direction of multiplet index!
%if isequal(para.sym,'SU2')
    B = getIdentityQS(phi(N),2,'-0');
    phi(N)=QSpace(contractQS(phi(N),[2],B,1,[1,3,2]));
    phi(N).info.itags = { strcat('B',num2str(N-1)),  strcat('B',num2str(N)) , strcat('s',num2str(N))};
%end


% we normalize the state
ovl = mpsOverlapQS(phi,phi);
IPhi.ovl = ovl;

if abs(ovl) <= 1e-16,
    printf('Info: The overlap <%s phi| %s phi> is zero',Fop,Fop);
    phi = QSpace;
    return;
end

phi(N) = phi(N) * (1/sqrt(ovl));


Phi = phi;
IPhi.block_idx = block_idx;
IPhi.state_idx = state_idx;

end

