function [A,Z,Info] = TrotterEvol(A,op,varargin)
%Carries out one (!) time step on a MPS A using a second order Trotter decomposition
%for a local NN Hamiltonian
%Input:
%- A: initial MPS to be time-evolved
%- op: structure containing local Trotter gates (op.TGates) for each bond
%- varargin:
%- varargin: options to be passed
%            'Nkeep': number of kept state in SVD
%            'Nt': number of time steps (default Nt=1)
%            'stol': minimum singular value kept in SVD
%            '-q': quiet-mode
%Output:
%- A: optimized MPS
%- Z: normalization weigth (imaginary time-evolution)
%- Info: structure containing information about truncation ect

%BB 2013, 11-2016; QSpace 3.0

global para

%input arguments
getopt('INIT',varargin);
qflag  = getopt('-q');
%check if time evolution is done in real time or imaginary time
if isreal(op.TGate) imagflag = 1; else imagflag=0; end

if ~isempty(find(strcmp('stol', varargin))) stol=getopt('stol',varargin); else stol=para.stol; end
if ~isempty(find(strcmp('Nkeep', varargin))) Nkeep=getopt('Nkeep',varargin); else Nkeep=para.Nkeep(1); end
if ~isempty(find(strcmp('Nt', varargin))) Nt=getopt('Nt',varargin); else Nt=1; end


%Combining two local indices to one enlarged index
if isfield(op,'FC') ID = QSpace(getIdentityQS(op.Z,1,op.Z,1)); elseif isfield(op,'S') ID = QSpace(getIdentityQS(op.S(1),1,op.S(1),1)); end

%%%%%% 2. order Trotter %%%%%%%
%first half-step (odd bonds)
for k=1:para.N-1
    if ~isequal(mod(k,2),0)
        if (length(A(k).Q) == 3 && length(A(k+1).Q) ==3)
            Q = contractQS(A(k),A(k+1)); % L s1 R s2
            Q = contractQS(Q,op.TGate2(k),[1,3,2,4]); %L s1 R s2
            %form bond matrix
            E1 = QSpace(getIdentityQS(Q,1,Q,2)); E1.info.itags{3} = 'L*';
            E2 = QSpace(getIdentityQS(Q,3,Q,4)); E2.info.itags{3} = 'R*';
            Q = contractQS(conj(E1),contractQS(Q,conj(E2)));
            %orthonormalization step
            opts = {'stol',stol,'Nkeep',Nkeep};
            [A(k+1), A(k), I1] = orthoQS(Q, 1, '>>', opts{:}); % A(k): sLR, A(k+1): LsR
            A(k).info.itags{1} =  ['B',num2str(k),'*'];
            A(k+1).info.itags{1} =  ['B',num2str(k)];
            A(k)=QSpace(contractQS(E1,A(k), [1 3 2]));
            A(k+1)=QSpace(contractQS(A(k+1),E2));
        elseif (length(A(k).Q)==4 && length(A(k+1).Q)==3)
            Q = contractQS(A(k),A(k+1)); % L s1 R s2
            Q = contractQS(Q,op.TGate2(k),[1,4,2,3,5]); %L s1 m R s2
            %form bond matrix
            E0 = getIdentityQS(Q,2,Q,3);
            E1 = getIdentityQS(Q,1,E0,3);
            E1 = QSpace(contractQS(E1,'2',E0,'3',[1,3,4,2])); E1.info.itags{4} = 'L*';
            E2 = QSpace(getIdentityQS(Q,4,Q,5)); E2.info.itags{3} = 'R*';
            Q = contractQS(conj(E1),contractQS(Q,conj(E2)));
            %orthonormalization step
            opts = {'stol',stol,'Nkeep',Nkeep};
            [A(k+1), A(k), I1] = orthoQS(Q, 1, '>>', opts{:}); % A(k): (Lsm) R, A(k+1): LsR
            A(k).info.itags{1} =  ['B',num2str(k),'*'];
            A(k+1).info.itags{1} =  ['B',num2str(k)];
            A(k)=QSpace(contractQS(E1,A(k), [1 4 2 3]));
            A(k+1)=QSpace(contractQS(A(k+1),E2));
        elseif (length(A(k).Q)==3 && length(A(k+1).Q)==4)
            Q = contractQS(A(k),A(k+1)); % L s1 R s2 m
            Q = contractQS(Q,op.TGate2(k),[1,4,3,5,2]); %L s1 R s2 m
            %form bond matrix
            E1 = QSpace(getIdentityQS(Q,1,Q,2)); E1.info.itags{3} = 'L*';
            E0 = getIdentityQS(Q,4,Q,5);
            E2 = getIdentityQS(Q,3,E0,3);
            E2 = QSpace(contractQS(E2,'2',E0,'3',[1,3,4,2])); E2.info.itags{4} = 'R*';
            Q = contractQS(conj(E1),contractQS(Q,conj(E2)));
            %orthonormalization step
            opts = {'stol',stol,'Nkeep',Nkeep};
            [A(k+1), A(k), I1] = orthoQS(Q, 1, '>>', opts{:}); % A(k): sLR, A(k+1): LsR
            A(k).info.itags{1} =  ['B',num2str(k),'*'];
            A(k+1).info.itags{1} =  ['B',num2str(k)];
            A(k)=QSpace(contractQS(E1,A(k), [1 3 2]));
            A(k+1)=QSpace(contractQS(A(k+1),E2));
        end
    else
        opts = {'stol',stol,'Nkeep',Nkeep,'nb',k,'-o'};
        [A(k),A(k+1),I1]=MPSOrthoQS(A(k),A(k+1),'>>',opts{:});
    end
end %for


%in case we do multiple time steps:
%group intermediate half steps
for n=2:Nt
    %intermediate full time step (even bonds)
    for k=para.N-1:-1:1
        if isequal(mod(k,2),0)
            if (length(A(k).Q) == 3 && length(A(k+1).Q) ==3)
                Q = contractQS(A(k),A(k+1)); % L s1 R s2
                Q = contractQS(Q,op.TGate(k),[1,3,2,4]); %L s1 R s2
                %form bond matrix
                E1 = QSpace(getIdentityQS(Q,1,Q,2)); E1.info.itags{3} = 'L*';
                E2 = QSpace(getIdentityQS(Q,3,Q,4)); E2.info.itags{3} = 'R*';
                Q = contractQS(conj(E1),contractQS(Q,conj(E2)));
                %orthonormalization step
                opts = {'stol',stol,'Nkeep',Nkeep};
                [A(k+1), A(k), I1] = orthoQS(Q, 1, '<<', opts{:}); % A(k): sLR, A(k+1): LsR
                A(k).info.itags{1} =  ['B',num2str(k)];
                A(k+1).info.itags{1} =  ['B',num2str(k),'*'];
                A(k)=QSpace(contractQS(E1,A(k), [1 3 2]));
                A(k+1)=QSpace(contractQS(A(k+1),E2));
            elseif (length(A(k).Q)==4 && length(A(k+1).Q)==3)
                Q = contractQS(A(k),A(k+1)); % L s1 R s2
                Q = contractQS(Q,op.TGate(k),[1,4,2,3,5]); %L s1 m R s2
                %form bond matrix
                E0 = getIdentityQS(Q,2,Q,3);
                E1 = getIdentityQS(Q,1,E0,3);
                E1 = QSpace(contractQS(E1,'2',E0,'3',[1,3,4,2])); E1.info.itags{4} = 'L*';
                E2 = QSpace(getIdentityQS(Q,4,Q,5)); E2.info.itags{3} = 'R*';
                Q = contractQS(conj(E1),contractQS(Q,conj(E2)));
                %orthonormalization step
                opts = {'stol',stol,'Nkeep',Nkeep};
                [A(k+1), A(k), I1] = orthoQS(Q, 1, '<<', opts{:}); % A(k): (Lsm) R, A(k+1): LsR
                A(k).info.itags{1} =  ['B',num2str(k)];
                A(k+1).info.itags{1} =  ['B',num2str(k),'*'];
                A(k)=QSpace(contractQS(E1,A(k), [1 4 2 3]));
                A(k+1)=QSpace(contractQS(A(k+1),E2));
            elseif (length(A(k).Q)==3 && length(A(k+1).Q)==4)
                Q = contractQS(A(k),A(k+1)); % L s1 R s2 m
                Q = contractQS(Q,op.TGate(k),[1,4,3,5,2]); %L s1 R s2 m
                %form bond matrix
                E1 = QSpace(getIdentityQS(Q,1,Q,2)); E1.info.itags{3} = 'L*';
                E0 = getIdentityQS(Q,4,Q,5);
                E2 = getIdentityQS(Q,3,E0,3);
                E2 = QSpace(contractQS(E2,'2',E0,'3',[1,3,4,2])); E2.info.itags{4} = 'R*';
                Q = contractQS(conj(E1),contractQS(Q,conj(E2)));
                %orthonormalization step
                opts = {'stol',stol,'Nkeep',Nkeep};
                [A(k+1), A(k), I1] = orthoQS(Q, 1, '<<', opts{:}); % A(k): sLR, A(k+1): LsR
                A(k).info.itags{1} =  ['B',num2str(k)];
                A(k+1).info.itags{1} =  ['B',num2str(k),'*'];
                A(k)=QSpace(contractQS(E1,A(k), [1 3 2]));
                A(k+1)=QSpace(contractQS(A(k+1),E2));
            end
        else
            opts = {'stol',stol,'Nkeep',Nkeep,'nb',k,'-o'};
            [A(k),A(k+1),I1]=MPSOrthoQS(A(k),A(k+1),'<<',opts{:});
        end
    end %for
    
    %intermediate full time step (odd bonds)
    for k=1:para.N-1
        if ~isequal(mod(k,2),0)
            if (length(A(k).Q) == 3 && length(A(k+1).Q) ==3)
                Q = contractQS(A(k),A(k+1)); % L s1 R s2
                Q = contractQS(Q,op.TGate(k),[1,3,2,4]); %L s1 R s2
                %form bond matrix
                E1 = QSpace(getIdentityQS(Q,1,Q,2)); E1.info.itags{3} = 'L*';
                E2 = QSpace(getIdentityQS(Q,3,Q,4)); E2.info.itags{3} = 'R*';
                Q = contractQS(conj(E1),contractQS(Q,conj(E2)));
                %orthonormalization step
                opts = {'stol',stol,'Nkeep',Nkeep};
                [A(k+1), A(k), I1] = orthoQS(Q, 1, '>>', opts{:}); % A(k): sLR, A(k+1): LsR
                A(k).info.itags{1} =  ['B',num2str(k),'*'];
                A(k+1).info.itags{1} =  ['B',num2str(k)];
                A(k)=QSpace(contractQS(E1,A(k), [1 3 2]));
                A(k+1)=QSpace(contractQS(A(k+1),E2));
            elseif (length(A(k).Q)==4 && length(A(k+1).Q)==3)
                Q = contractQS(A(k),A(k+1)); % L s1 R s2
                Q = contractQS(Q,op.TGate(k),[1,4,2,3,5]); %L s1 m R s2
                %form bond matrix
                E0 = getIdentityQS(Q,2,Q,3);
                E1 = getIdentityQS(Q,1,E0,3);
                E1 = QSpace(contractQS(E1,'2',E0,'3',[1,3,4,2])); E1.info.itags{4} = 'L*';
                E2 = QSpace(getIdentityQS(Q,4,Q,5)); E2.info.itags{3} = 'R*';
                Q = contractQS(conj(E1),contractQS(Q,conj(E2)));
                %orthonormalization step
                opts = {'stol',stol,'Nkeep',Nkeep};
                [A(k+1), A(k), I1] = orthoQS(Q, 1, '>>', opts{:}); % A(k): (Lsm) R, A(k+1): LsR
                A(k).info.itags{1} =  ['B',num2str(k),'*'];
                A(k+1).info.itags{1} =  ['B',num2str(k)];
                A(k)=QSpace(contractQS(E1,A(k), [1 4 2 3]));
                A(k+1)=QSpace(contractQS(A(k+1),E2));
            elseif (length(A(k).Q)==3 && length(A(k+1).Q)==4)
                Q = contractQS(A(k),A(k+1)); % L s1 R s2 m
                Q = contractQS(Q,op.TGate(k),[1,4,3,5,2]); %L s1 R s2 m
                %form bond matrix
                E1 = QSpace(getIdentityQS(Q,1,Q,2)); E1.info.itags{3} = 'L*';
                E0 = getIdentityQS(Q,4,Q,5);
                E2 = getIdentityQS(Q,3,E0,3);
                E2 = QSpace(contractQS(E2,'2',E0,'3',[1,3,4,2])); E2.info.itags{4} = 'R*';
                Q = contractQS(conj(E1),contractQS(Q,conj(E2)));
                %orthonormalization step
                opts = {'stol',stol,'Nkeep',Nkeep};
                [A(k+1), A(k), I1] = orthoQS(Q, 1, '>>', opts{:}); % A(k): sLR, A(k+1): LsR
                A(k).info.itags{1} =  ['B',num2str(k),'*'];
                A(k+1).info.itags{1} =  ['B',num2str(k)];
                A(k)=QSpace(contractQS(E1,A(k), [1 3 2]));
                A(k+1)=QSpace(contractQS(A(k+1),E2));
            end
        else
            opts = {'stol',stol,'Nkeep',Nkeep,'nb',k,'-o'};
            [A(k),A(k+1),I1]=MPSOrthoQS(A(k),A(k+1),'>>',opts{:});
        end
    end %for
    
end



Info.svdtr=[]; Info.Nkeep=[];
%backward sweep - evolve even bonds!
for k=para.N-1:-1:1
    if isequal(mod(k,2),0)
        if (length(A(k).Q) == 3 && length(A(k+1).Q) ==3)
            Q = contractQS(A(k),A(k+1)); % L s1 R s2
            Q = contractQS(Q,op.TGate(k),[1,3,2,4]); %L s1 R s2
            %form bond matrix
            E1 = QSpace(getIdentityQS(Q,1,Q,2)); E1.info.itags{3} = 'L*';
            E2 = QSpace(getIdentityQS(Q,3,Q,4)); E2.info.itags{3} = 'R*';
            Q = contractQS(conj(E1),contractQS(Q,conj(E2)));
            %orthonormalization step
            opts = {'stol',stol,'Nkeep',Nkeep};
            [A(k+1), A(k), I1] = orthoQS(Q, 1, '<<', opts{:}); % A(k): sLR, A(k+1): LsR
            A(k).info.itags{1} =  ['B',num2str(k)];
            A(k+1).info.itags{1} =  ['B',num2str(k),'*'];
            A(k)=QSpace(contractQS(E1,A(k), [1 3 2]));
            A(k+1)=QSpace(contractQS(A(k+1),E2));
        elseif (length(A(k).Q)==4 && length(A(k+1).Q)==3)
            Q = contractQS(A(k),A(k+1)); % L s1 R s2
            Q = contractQS(Q,op.TGate(k),[1,4,2,3,5]); %L s1 m R s2
            %form bond matrix
            E0 = getIdentityQS(Q,2,Q,3);
            E1 = getIdentityQS(Q,1,E0,3);
            E1 = QSpace(contractQS(E1,'2',E0,'3',[1,3,4,2])); E1.info.itags{4} = 'L*';
            E2 = QSpace(getIdentityQS(Q,4,Q,5)); E2.info.itags{3} = 'R*';
            Q = contractQS(conj(E1),contractQS(Q,conj(E2)));
            %orthonormalization step
            opts = {'stol',stol,'Nkeep',Nkeep};
            [A(k+1), A(k), I1] = orthoQS(Q, 1, '<<', opts{:}); % A(k): (Lsm) R, A(k+1): LsR
            A(k).info.itags{1} =  ['B',num2str(k)];
            A(k+1).info.itags{1} =  ['B',num2str(k),'*'];
            A(k)=QSpace(contractQS(E1,A(k), [1 4 2 3]));
            A(k+1)=QSpace(contractQS(A(k+1),E2));
        elseif (length(A(k).Q)==3 && length(A(k+1).Q)==4)
            Q = contractQS(A(k),A(k+1)); % L s1 R s2 m
            Q = contractQS(Q,op.TGate(k),[1,4,3,5,2]); %L s1 R s2 m
            %form bond matrix
            E1 = QSpace(getIdentityQS(Q,1,Q,2)); E1.info.itags{3} = 'L*';
            E0 = getIdentityQS(Q,4,Q,5);
            E2 = getIdentityQS(Q,3,E0,3);
            E2 = QSpace(contractQS(E2,'2',E0,'3',[1,3,4,2])); E2.info.itags{4} = 'R*';
            Q = contractQS(conj(E1),contractQS(Q,conj(E2)));
            %orthonormalization step
            opts = {'stol',stol,'Nkeep',Nkeep};
            [A(k+1), A(k), I1] = orthoQS(Q, 1, '<<', opts{:}); % A(k): sLR, A(k+1): LsR
            A(k).info.itags{1} =  ['B',num2str(k)];
            A(k+1).info.itags{1} =  ['B',num2str(k),'*'];
            A(k)=QSpace(contractQS(E1,A(k), [1 3 2]));
            A(k+1)=QSpace(contractQS(A(k+1),E2));
        end
    else
        opts = {'stol',stol,'Nkeep',Nkeep,'nb',k,'-o'};
        [A(k),A(k+1),I1]=MPSOrthoQS(A(k),A(k+1),'<<',opts{:});
    end
end %for


%last half-step -odd bonds
for k=1:para.N-1
    if ~isequal(mod(k,2),0)
        if (length(A(k).Q) == 3 && length(A(k+1).Q) ==3)
            Q = contractQS(A(k),A(k+1)); % L s1 R s2
            Q = contractQS(Q,op.TGate2(k),[1,3,2,4]); %L s1 R s2
            %form bond matrix
            E1 = QSpace(getIdentityQS(Q,1,Q,2)); E1.info.itags{3} = 'L*';
            E2 = QSpace(getIdentityQS(Q,3,Q,4)); E2.info.itags{3} = 'R*';
            Q = contractQS(conj(E1),contractQS(Q,conj(E2)));
            %orthonormalization step
            opts = {'stol',stol,'Nkeep',Nkeep};
            [A(k+1), A(k), I1] = orthoQS(Q, 1, '>>', opts{:}); % A(k): sLR, A(k+1): LsR
            A(k).info.itags{1} =  ['B',num2str(k),'*'];
            A(k+1).info.itags{1} =  ['B',num2str(k)];
            A(k)=QSpace(contractQS(E1,A(k), [1 3 2]));
            A(k+1)=QSpace(contractQS(A(k+1),E2));
            Info.svdtr=[Info.svdtr,I1.svd2tr];
            Info.Nkeep=[Info.Nkeep,I1.Nkeep];
        elseif (length(A(k).Q)==4 && length(A(k+1).Q)==3)
            Q = contractQS(A(k),A(k+1)); % L s1 R s2
            Q = contractQS(Q,op.TGate2(k),[1,4,2,3,5]); %L s1 m R s2
            %form bond matrix
            E0 = getIdentityQS(Q,2,Q,3);
            E1 = getIdentityQS(Q,1,E0,3);
            E1 = QSpace(contractQS(E1,'2',E0,'3',[1,3,4,2])); E1.info.itags{4} = 'L*';
            E2 = QSpace(getIdentityQS(Q,4,Q,5)); E2.info.itags{3} = 'R*';
            Q = contractQS(conj(E1),contractQS(Q,conj(E2)));
            %orthonormalization step
            opts = {'stol',stol,'Nkeep',Nkeep};
            [A(k+1), A(k), I1] = orthoQS(Q, 1, '>>', opts{:}); % A(k): (Lsm) R, A(k+1): LsR
            A(k).info.itags{1} =  ['B',num2str(k),'*'];
            A(k+1).info.itags{1} =  ['B',num2str(k)];
            A(k)=QSpace(contractQS(E1,A(k), [1 4 2 3]));
            A(k+1)=QSpace(contractQS(A(k+1),E2));
            Info.svdtr=[Info.svdtr,I1.svd2tr];
            Info.Nkeep=[Info.Nkeep,I1.Nkeep];
        elseif (length(A(k).Q)==3 && length(A(k+1).Q)==4)
            Q = contractQS(A(k),A(k+1)); % L s1 R s2 m
            Q = contractQS(Q,op.TGate2(k),[1,4,3,5,2]); %L s1 R s2 m
            %form bond matrix
            E1 = QSpace(getIdentityQS(Q,1,Q,2)); E1.info.itags{3} = 'L*';
            E0 = getIdentityQS(Q,4,Q,5);
            E2 = getIdentityQS(Q,3,E0,3);
            E2 = QSpace(contractQS(E2,'2',E0,'3',[1,3,4,2])); E2.info.itags{4} = 'R*';
            Q = contractQS(conj(E1),contractQS(Q,conj(E2)));
            %orthonormalization step
            opts = {'stol',stol,'Nkeep',Nkeep};
            [A(k+1), A(k), I1] = orthoQS(Q, 1, '>>', opts{:}); % A(k): sLR, A(k+1): LsR
            A(k).info.itags{1} =  ['B',num2str(k),'*'];
            A(k+1).info.itags{1} =  ['B',num2str(k)];
            A(k)=QSpace(contractQS(E1,A(k), [1 3 2]));
            A(k+1)=QSpace(contractQS(A(k+1),E2));
            Info.svdtr=[Info.svdtr,I1.svd2tr];
            Info.Nkeep=[Info.Nkeep,I1.Nkeep];
        end
    else
        opts = {'stol',stol,'Nkeep',Nkeep,'nb',k,'-o'};
        [A(k),A(k+1),I1]=MPSOrthoQS(A(k),A(k+1),'>>',opts{:});
        Info.svdtr=[Info.svdtr,I1.svd2tr];
        Info.Nkeep=[Info.Nkeep,I1.Nkeep];
    end
end %for



for k=para.N-1:-1:1
    opts = {'stol',stol,'Nkeep',Nkeep,'nb',k,'-o'};
    [A(k),A(k+1),I1]=MPSOrthoQS(A(k),A(k+1),'<<',opts{:});
end

%only normalize in case of imaginary time evolution
if imagflag
    Z = normQS(A(1));
    A(1) = A(1)*(1/Z);
else
    Z=1;
end

end %function
