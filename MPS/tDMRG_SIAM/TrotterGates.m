function op = TrotterGates(op,tau)
%generates two-site Trotter gates for NN-Fermionic Hamiltonian
%Input:
%- op: structure containing local operators
%- tau: time step (real or imaginary!) 
%Output:
%- op: operator structure now including Trotter gates

%BB 2013, 11-2016; QSpace 3.0

global para


switch para.model
    case 'FermionS_NN'
        op.TGate=QSpace(para.N-1,1);
        
        for k=1:para.N-1
            %write bond Hamiltonian h in bond basis
            ID = QSpace(getIdentityQS(op.h1{k,1},1,op.h1{k+1,1},1));
            ID.info.itags{3} = ['s',num2str(k),'s',num2str(k+1),'*'];
            Q=QSpace;
            h=QSpace;
                     
            for l=1:length(op.h2r(k,:))
                %COUPLING TERMS
                Q=contractQS(contractQS(op.h2r{k,l},ID),op.h2l{k+1,l},[1,3,2]);
                Q=QSpace(contractQS(ID,'12*',Q,'12'));
                h=h+Q;
            end
            
            %onsite terms (including local interaction!)
            Q=contractQS(op.h1{k,1},'2',ID,'1');
            Q=QSpace(contractQS(ID,'12*',Q,'12'));
            h=h+Q;
      
            if isequal(k,para.N-1)
                     %onsite terms (including local interaction!)
                     Q=contractQS(op.h1{k+1,1},'2',ID,'2',[2,1,3]);
                     Q=QSpace(contractQS(ID,'12*',Q,'12'));
                    h=h+Q;
            end
            
            %add missing QSpaces
            Q=QSpace(contractQS(1e-25*ID,'12*',ID,'12'));
            h=h+Q;
            
            %exponentiate
            for l=1:length(h.data)
                h.data{l} = expm(-tau.*h.data{l});
            end
            h = QSpace(h);
%             h = contractQS(conj(ID),3,h,2);
%             h = contractQS(ID,3,h,3);

h= contractQS(ID,'3',h,'1');
h = QSpace(contractQS(h,3,ID,'3*'));            


            op.TGate(k) = QSpace(h);
            
            
            
        end
     
        
        
end %switch model


end% function

