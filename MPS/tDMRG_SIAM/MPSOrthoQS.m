function [A1,A2,I] = MPSOrthoQS(A1,A2,lrdir,varargin)
%orthogonalization of MPS with non-abelian symmetries
%Input:
%- A1,A2: describe neighbouring A matrices of the MPS
%- lrdir: sets direction for canonical update ('>>' for right-move, '<<' for left-move)
%- varargin: options to be passed
%	     'Nkeep': number of kept state in SVD 
%	     'stol': minimum singular value kept in SVD 
%	     'nb': label (number) of bond in MPS to be updated 
%	     '-o': orthogonalization flag (if set, plain orthonormalization is performed; otherwise, full two-site tensor is formed) 

%BB 01-2015, 10-2015, 11-2016; QSpace 3.0

global para

%input options (see above)
getopt('INIT',varargin);
oflag  = getopt('-o'); 
if ~isempty(find(strcmp('stol', varargin))) stol=getopt('stol',varargin); else stol=para.stol; end
if ~isempty(find(strcmp('Nkeep', varargin))) Nkeep=getopt('Nkeep',varargin); else Nkeep=para.Nkeep(1); end
if ~isempty(find(strcmp('nb', varargin))) nb = getopt('nb',varargin); end

%set options for SVD (see above)
opts = {'Nkeep',Nkeep,'stol',stol};

%dimension check
d1 = length(A1.Q); d2 = length(A2.Q);

%standard: both tensors have three indices
if (d1==3 && d2==3)
    if oflag
        if isequal(lrdir,'>>')
            [U,A1,I] =  orthoQS(A1,[1,3],'>>',opts{:});
            A1 = QSpace(A1);
            A2 = QSpace(contractQS(U,1,A2,1));
            %set new itags
            if exist('nb','var')
                A1.info.itags = { ['B',num2str(nb-1)],  ['B',num2str(nb),'*'] , ['s',num2str(nb)]};
                A2.info.itags = { ['B',num2str(nb)],  ['B',num2str(nb+1)] , ['s',num2str(nb+1)]};
            end
        elseif isequal(lrdir,'<<')
            [A2,U,I] =  orthoQS(A2,1,'<<',opts{:});
            A2 = QSpace(A2);
            A1 = QSpace(contractQS(A1,2,U,2,[1,3,2]));
            %set new itags
            if exist('nb','var')
                A1.info.itags = { ['B',num2str(nb-1)],  ['B',num2str(nb)] , ['s',num2str(nb)]};
                A2.info.itags = { ['B',num2str(nb),'*'],  ['B',num2str(nb+1)] , ['s',num2str(nb+1)]};
            end
        end%lrdir
    else
        A=contractQS(A1,2,A2,1);% l s1 r s2
        %isometries
        isom1 = QSpace(getIdentityQS(A,1,A,2));
        isom2 = QSpace(getIdentityQS(A,3,A,4));
        isom1.info.itags{3} = 'cl*';
        isom2.info.itags{3} = 'cr*';
        %projecting 2-site block in effective basis
        A=contractQS(isom1,'1 2;*',A,'1 2');
        A=QSpace(contractQS(A,'2 3',isom2,'1 2;*')); %(l s1) (r s2)
        if isequal(lrdir,'>>')
            [A2,A1,I]=orthoQS(A,1,'>>',opts{:});
            A1 = QSpace(contractQS(isom1,3,A1,2,[1,3,2]));
            A2 = QSpace(contractQS(A2,2,isom2,[3]));
            %set new itags
            if exist('nb','var')
                A1.info.itags = { ['B',num2str(nb-1)],  ['B',num2str(nb),'*'] , ['s',num2str(nb)]};
                A2.info.itags = { ['B',num2str(nb)],  ['B',num2str(nb+1)] , ['s',num2str(nb+1)]};
            end
        elseif isequal(lrdir,'<<')
            [A2,A1,I]=orthoQS(A,1,'<<',opts{:});
            A1 = QSpace(contractQS(isom1,3,A1,2,[1,3,2]));
            A2 = QSpace(contractQS(A2,2,isom2,[3]));
            %set new itags
            if exist('nb','var')
                A1.info.itags = { ['B',num2str(nb-1)],  ['B',num2str(nb)] , ['s',num2str(nb)]};
                A2.info.itags = { ['B',num2str(nb),'*'],  ['B',num2str(nb+1)] , ['s',num2str(nb+1)]};
            end
        end%lrdir
    end %oflag
    
%non-standard: one tensor has additional index due to presence of
%excitation (e.g. A1" = c'A1)
elseif (d1==4 && d2==3)
    itag_m = A1.info.itags{4};
    if oflag
        if isequal(lrdir,'>>')
            [U,A1,I] =  orthoQS(A1,[1,3,4],'>>',opts{:});
            A1 = QSpace(A1);
            A2 = QSpace(contractQS(U,1,A2,1));
            %set new itags
            if exist('nb','var')
                A1.info.itags = { ['B',num2str(nb-1)],  ['B',num2str(nb),'*'] , ['s',num2str(nb)], itag_m};
                A2.info.itags = { ['B',num2str(nb)],  ['B',num2str(nb+1)] , ['s',num2str(nb+1)]};
            end
        elseif isequal(lrdir,'<<')
            [A2,U,I] =  orthoQS(A2,1,'<<',opts{:});
            A2 = QSpace(A2);
            A1 = QSpace(contractQS(A1,2,U,2,[1,4,2,3]));
            %set new itags
            if exist('nb','var')
                A1.info.itags = { ['B',num2str(nb-1)],  ['B',num2str(nb)] , ['s',num2str(nb)],itag_m};
                A2.info.itags = { ['B',num2str(nb),'*'],  ['B',num2str(nb+1)] , ['s',num2str(nb+1)]};
            end
        end%lrdir
    else
        A=contractQS(A1,2,A2,1);% l s1 m r s2
        %isometries
        isom0 = getIdentityQS(A,2,A,3);
        isom1 = getIdentityQS(A,1,isom0,3);
        isom1 = QSpace(contractQS(isom1,2,isom0,3,[1,3,4,2]));
        isom2 = QSpace(getIdentityQS(A,4,A,5));
        isom1.info.itags{4} = 'cl*';
        isom2.info.itags{3} = 'cr*';
        %projecting 2-site block in effective basis
        A=contractQS(isom1,'1 2 3;*',A,'1 2 3');
        A=QSpace(contractQS(A,'2 3',isom2,'1 2;*')); %(l s1 m) (r s2)
        if isequal(lrdir,'>>')
            [A2,A1,I]=orthoQS(A,1,'>>',opts{:});
            A1 = QSpace(contractQS(isom1,'4',A1,'2',[1,4,2,3]));
            A2 = QSpace(contractQS(A2,'2',isom2,'3'));
            %set new itags
            if exist('nb','var')
                A1.info.itags = { ['B',num2str(nb-1)],  ['B',num2str(nb),'*'] , ['s',num2str(nb)],itag_m};
                A2.info.itags = { ['B',num2str(nb)],  ['B',num2str(nb+1)] , ['s',num2str(nb+1)]};
            end
        elseif isequal(lrdir,'<<')
            [A2,A1,I]=orthoQS(A,1,'<<',opts{:});
            A1 = QSpace(contractQS(isom1,'4',A1,'2',[1,4,2,3]));
            A2 = QSpace(contractQS(A2,'2',isom2,'3'));
            %set new itags
            if exist('nb','var')
                A1.info.itags = { ['B',num2str(nb-1)],  ['B',num2str(nb)] , ['s',num2str(nb)],itag_m};
                A2.info.itags = { ['B',num2str(nb),'*'],  ['B',num2str(nb+1)] , ['s',num2str(nb+1)]};
            end
        end%lrdir
    end %oflag
    
elseif (d1==3 && d2==4)
    itag_m = A2.info.itags{4};
    if oflag
        if isequal(lrdir,'>>')
            [U,A1,I] =  orthoQS(A1,[1,3],'>>',opts{:});
            A1 = QSpace(A1);
            A2 = QSpace(contractQS(U,1,A2,1));
            %set new itags
            if exist('nb','var')
                A1.info.itags = { ['B',num2str(nb-1)],  ['B',num2str(nb),'*'] , ['s',num2str(nb)]};
                A2.info.itags = { ['B',num2str(nb)],  ['B',num2str(nb+1)] , ['s',num2str(nb+1)],itag_m};
            end
        elseif isequal(lrdir,'<<')
            [A2,U,I] =  orthoQS(A2,1,'<<',opts{:});
            A2 = QSpace(A2);
            A1 = QSpace(contractQS(A1,2,U,2,[1,3,2]));
            %set new itags
            if exist('nb','var')
                A1.info.itags = { ['B',num2str(nb-1)],  ['B',num2str(nb)] , ['s',num2str(nb)]};
                A2.info.itags = { ['B',num2str(nb),'*'],  ['B',num2str(nb+1)] , ['s',num2str(nb+1)],itag_m};
            end
        end%lrdir
    else
        A=contractQS(A1,2,A2,1);% l s1 r s2
        %isometries
        isom1 = QSpace(getIdentityQS(A,'1',A,'2'));
        isom0 = getIdentityQS(A,'4',A,'5');
        isom2 = getIdentityQS(A,'3',isom0,'3');
        isom2 = QSpace(contractQS(isom2,'2',isom0,'3',[1,3,4,2]));
        isom1.info.itags{3} = 'cl*';
        isom2.info.itags{4} = 'cr*';
        %projecting 2-site block in effective basis
        A=contractQS(isom1,'1 2;*',A,'1 2');
        A=QSpace(contractQS(A,'2 3 4',isom2,'1 2 3;*')); %(l s1) (r s2 m)
        if isequal(lrdir,'>>')
            [A2,A1,I]=orthoQS(A,1,'>>',opts{:});
            A1 = QSpace(contractQS(isom1,'3',A1,'2',[1,3,2]));
            A2 = QSpace(contractQS(A2,'2',isom2,'4'));
            %set new itags
            if exist('nb','var')
                A1.info.itags = { ['B',num2str(nb-1)],  ['B',num2str(nb),'*'] , ['s',num2str(nb)]};
                A2.info.itags = { ['B',num2str(nb)],  ['B',num2str(nb+1)] , ['s',num2str(nb+1)],itag_m};
            end
        elseif isequal(lrdir,'<<')
            [A2,A1,I]=orthoQS(A,1,'<<',opts{:});
            A1 = QSpace(contractQS(isom1,'3',A1,'2',[1,3,2]));
            A2 = QSpace(contractQS(A2,'2',isom2,'4'));
            %set new itags
            if exist('nb','var')
                A1.info.itags = { ['B',num2str(nb-1)],  ['B',num2str(nb)] , ['s',num2str(nb)]};
                A2.info.itags = { ['B',num2str(nb),'*'],  ['B',num2str(nb+1)] , ['s',num2str(nb+1)],itag_m};
            end
        end%lrdir
    end %oflag
       
else
    error('ERR: A1 or A2 have not the correct number of indices!');
end

end%function
