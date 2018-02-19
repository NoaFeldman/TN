function op=OpUpdate(A,op,k,lrdir)
% updates global left and right Hamiltonian parts (HLR) as well as local coupling terms (CLR) after othogonalization step to neighbouring site
%Input:
%- A: describes A matrices of the MPS at site k (assumed to be right- or left-canonical!)
%- lrdir: sets direction for canonical update ('>>' for right-move, '<<' for left-move)
%- op: structure containing local Hamiltonian terms as well as global left and right Hamiltonian parts
%Output:
%- op: updated operator structure containing local Hamiltonian terms as well as global left and right Hamiltonian parts

% BB 01-2013, 11-2016; QSpace 3.0


global para 

%number of coupling terms
m=size(op.h2r,2);


%%%% UPDATE WHILE SWEEPING FROM R TO L %%%%
if isequal(lrdir,'<<')
    Q=QSpace;
    %setting up next neighbour coupling
    for i=1:m
        Q=contractQS(A,op.h2l{k,i});
        op.CLR{k,i}=QSpace(contractQS(A,'23*',Q,'23'));
    end
    if ~isequal(k,para.N)    
        %HLR propagate to next site
        Q=contractQS(A,op.HLR(k+1)); % LsR
        op.HLR(k)=QSpace(contractQS(A,'23*',Q,'32'));
        %HLR: adding coupling terms
        for i=1:m
            Q=contractQS(contractQS(A,op.h2r{k,i}),op.CLR{k+1,i});
            Q=QSpace(contractQS(A,'23*',Q,'32'));
            op.HLR(k)=op.HLR(k)+Q;
        end
    else
       op.HLR(para.N)=QSpace;   
    end
    %HLR: adding local term
    if ~isempty(op.h1{k}) 
        Q=contractQS(A,op.h1{k});
        Q=contractQS(A,'23*',Q,'23');
        op.HLR(k)=op.HLR(k)+QSpace(Q);
    end


%%%% UPDATE WHILE SWEEPING FROM L TO R %%%%
elseif isequal(lrdir,'>>')
    Q=QSpace;
    %setting up next neighbour coupling
    for i=1:m
        Q=contractQS(A,op.h2r{k,i});
        op.CLR{k,i}=QSpace(contractQS(A,'13*',Q,'13'));
    end
    if isequal(k,1)
        %for site 1 the setup is the same to H0
        op.HLR(1)=QSpace;
        if ~isempty(op.h1{k})
            Q=contractQS(A,op.h1{k});
            Q=contractQS(A,'13*',Q,'13');
            op.HLR(k)=QSpace(Q);
        end
    else
        %HLR propagate to next site
        Q=contractQS(op.HLR(k-1),A); % LsR
        op.HLR(k)=QSpace(contractQS(A,'13*',Q,'13'));
        %HLR: adding coupling terms
        for i=1:m
           Q=contractQS(contractQS(op.CLR{k-1,i},A),op.h2l{k,i});
           Q=QSpace(contractQS(A,'13*',Q,'13'));
           op.HLR(k)=op.HLR(k)+Q;
        end
         %HLR: adding local term 
         if ~isempty(op.h1{k})
            Q=contractQS(A,op.h1{k});
            Q=contractQS(A,'13*',Q,'13');
            op.HLR(k)=op.HLR(k)+QSpace(Q);
         end
     end   
 
else fprint('\n Invalid lrdir %d\t',lrdir);
end
end
