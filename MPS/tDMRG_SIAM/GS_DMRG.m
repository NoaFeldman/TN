function [A,op,results]=GS_DMRG(A,op,results,varargin)
%standard DMRG optimization to determine ground state of system with nearest-neighbour interactions only in QSpace
%Input:
%- A: initial MPS to be optimized (vector of QSpaces)
%- op: structure containing local Hamiltonian terms as well as global left and right Hamiltonian parts
%- results: structure to save results such as ground-state energy
%- varargin: 
%- varargin: options to be passed
%            'Nkeep': number of kept state in SVD
%            'stol': minimum singular value kept in SVD
%            '-q': quiet-mode - information of Lanczos optimization not printed
%Output:
%- A: optimized MPS

%BB 11-2012, 11-2016; QSpace 3.0

global para


%input options (see above)
getopt('INIT',varargin);
qflag  = getopt('-q');
if ~isempty(find(strcmp('stol', varargin))) stol=getopt('stol',varargin); else stol=para.stol; end
if ~isempty(find(strcmp('Nkeep', varargin))) Nkeep=getopt('Nkeep',varargin); else Nkeep=para.Nkeep(1); end

%set options for SVD 
if qflag
    opts = {'-q'};
else
    opts = {};
end
opts = {opts{:},'Nkeep',para.Nkeep,'stol',para.stol};

% check which update procedure to use
switch para.update
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    case '2bond' %2-site bond update -> necessary when exploiting any kind of abelian or non-abelian symmetry
        results.SVDfull=cell(1,1);
        
        while para.sw <= para.swmax
            if isequal(para.sw,1)
                results.Eerror=1;
            end
            
            %%%%%%%%%%%%%%%
            %forward sweep%
            %%%%%%%%%%%%%%%
            fprintf('\v\n DMRG forward sweep %g \n', para.sw);
            fprintf(' ');
            for k=1:(para.N-1)
                fprintf('%g-', k);
                %optimization step
                [A(k),A(k+1),E,Isvd]=EigsQS_2bond(A(k),A(k+1),op,k,'>>',opts{:}); %includes orthonomalization step
                %operator update
                op=OpUpdate(A(k),op,k,'>>');
                results.E(k,para.sw)=E;
                results.SVDtr(k,para.sw)=Isvd.svd2tr;
                results.Nkeep_DMRG(k,para.sw) = Isvd.Nkeep;
                results.SVDfull{end+1,para.sw} = Isvd.svd;
            end
            
            %calculate energy error and output sweep info
            fprintf('\v \n Energy = %.10g\t', results.E(para.N-1,para.sw));
            results.Eerror(para.sw,1)=std(results.E(1:para.N-1,para.sw))/abs(mean(results.E(1:para.N-1,para.sw)));
            fprintf('\n Energy error =  %.16g\t',results.Eerror(para.sw,1));
            fprintf('\n Bond     Nkeep    SVDtr \n');
            for k=1:para.N-1
                fprintf(' %3.0f       %4.0f     %.5g \n', k, results.Nkeep_DMRG(k,para.sw), results.SVDtr(k,para.sw));
            end
            
            %%%%%%%%%%%%%%%%
            %backward sweep%
            %%%%%%%%%%%%%%%%
            fprintf('\v\n DMRG backward sweep %g \n', para.sw);
            fprintf(' ');
            for k=para.N-1:-1:1
                fprintf('%g-', k);
                %optimization step
                [A(k),A(k+1),E,Isvd]=EigsQS_2bond(A(k),A(k+1),op,k,'<<',opts{:}); %includes orthonomalization step
                %operator update
                op=OpUpdate(A(k+1),op,k+1,'<<');
                results.E(para.N-1+(para.N-(k)),para.sw)=E;
                results.SVDtr(para.N-1+(para.N-(k)),para.sw)=Isvd.svd2tr;
                results.Nkeep_DMRG(para.N-1+(para.N-(k)),para.sw) = Isvd.Nkeep;
                results.SVDfull{end+1,para.sw} = Isvd.svd;
            end
            
            %calculate energy error and output sweep info
            fprintf('\v \n Energy = %.10g\t', results.E(end,para.sw));
            results.Eerror(para.sw,2)=std(results.E(para.N:end,para.sw))/abs(mean(results.E(para.N:end,para.sw)));
            fprintf('\n Energy error =  %.16g\t',results.Eerror(para.sw,2));
            fprintf('\n Bond     Nkeep     SVDtr \n');
            for k=1:para.N-1
                fprintf(' %3.0f       %4.0f     %.5g \n', k, results.Nkeep_DMRG(para.N-1+(para.N-(k)),para.sw), results.SVDtr(para.N-1+(para.N-(k)),para.sw));
            end
            
            %convergence criterium
            if results.Eerror(para.sw,2)<para.precision
                break;
            end
            
%             if mod(para.sw,10)==0
%                 save(para.filename,'para','A','results','op');
%             end
            para.sw=para.sw+1;
            fprintf('\n %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% \n ');
        end %sweep while loop
end %switch update

end %function







