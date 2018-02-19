function op = GenHamiltonian(op,varargin)
%generates Hamiltonian for iterative DMRG optimization
%Input:
%- op: structure containing local site operators
%- varargin: currently not used; intended for rescaling ect
%Output:
%- op: updated structure additionally containing on-site and nearest-neighbour terms of Hamiltonian for specified model

%BB 01-2013, 11-2016; QSpace 3.0


global para

%input arguments
getopt('INIT',varargin);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%operator setup:
%op.h1{index_site}: local on-site terms
%op.h2l{index_site,index_cp}: coupling terms to the left (site k-1)
%op.h2r{index_site,index_cp}: coupling terms to the right (site k+1)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
switch para.model
    case 'FermionS_NN'
        %check that all parameters are set; use default values otherwise
        if ~isfield(para,'hz') para.hz=zeros(para.N,1); end
        if ~isfield(para,'V') para.V=zeros(para.N,1); end
        if ~isfield(para,'U') para.U=zeros(para.N,1); end
        if ~isfield(para,'t') para.t=ones(para.N,1); end
        %setup operators on different site
        for k=1:para.N
            if (isequal(para.sym,'Acharge,Aspin') | isequal(para.sym,'Pcharge,Aspin') | isequal(para.sym,'Aspin'))
                op.h1{k,1}=(-para.mu(k)+para.V(k)-para.hz(k)/2)*(op.N(1))+(-para.mu(k)+para.V(k)+para.hz(k)/2)*op.N(2)+para.U(k)*op.N(1)*op.N(2);
            elseif (isequal(para.sym,'Acharge,SU2spin') | isequal(para.sym,'SU2spin') | isequal(para.sym,'Pcharge,SU2spin'))
                op.h1{k,1}=(-para.mu(k)+para.V(k))*(op.N(1))+para.U(k)*(op.N(1)-op.ID)*op.N(1);
            elseif (isequal(para.sym,'SU2charge,Aspin') | isequal(para.sym,'SU2charge'))
                op.h1{k,1}=(-para.hz(k)/2)*(op.N(1))+(para.hz(k)/2)*op.N(2)+para.U(k)*op.N(1)*op.N(2);
            elseif isequal(para.sym,'SU2charge,SU2spin')
                op.h1{k,1}=+para.U(k)*(op.N(1)-op.ID)*op.N(1);
            end
            op.h1{k,1} = op.h1{k,1} + 1e-20*op.ID;
            if k==1
                for l = 1:length(op.FC)
                    op.h2l{k,2*l-1}=QSpace;      op.h2r{k,2*l-1}=-para.t(k)*op.FC(l);
                    op.h2l{k,2*l}=QSpace;        op.h2r{k,2*l}=-para.t(k)*op.FC(l)';
                end
            elseif k==para.N
                for l = 1:length(op.FC)
                    op.h2l{k,2*l-1}=(op.Z*op.FC(l))';  op.h2r{k,2*l-1}=QSpace;
                    op.h2l{k,2*l}=op.Z*op.FC(l);       op.h2r{k,2*l}=QSpace;
                end
            else
                for l = 1:length(op.FC)
                    op.h2l{k,2*l-1}=(op.Z*op.FC(l))';  op.h2r{k,2*l-1,1}=-para.t(k)*op.FC(l);
                    op.h2l{k,2*l}=op.Z*op.FC(l);       op.h2r{k,2*l}=-para.t(k)*op.FC(l)';
                end
            end % if k
            
            %label itags
            op.h1{k,1}.info.itags = {['s',num2str(k)],  ['s',num2str(k),'*']};
            for l = 1:length(op.h2l(k,:))
                op.h2l{k,l}.info.itags{1} = ['s',num2str(k)];  op.h2l{k,l}.info.itags{2} = ['s',num2str(k),'*'];
                op.h2r{k,l}.info.itags{1} = ['s',num2str(k)];  op.h2r{k,l}.info.itags{2} = ['s',num2str(k),'*'];
            end
        end %for k
end %switch model
       
end %function
