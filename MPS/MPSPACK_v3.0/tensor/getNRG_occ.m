function [nn,Iocc]=getNRG_occ(NN,varargin)
% function [nn,Iocc]=getNRG_occ(NN [,nrg,opts])
%
%   NN    set of local occupation operators to be used
%   nrg   NRG tag ($LMA/NRG/NRG)
%
% in addition, occupations regarding the Wilson sites are calculated
% with the data returned as nn(:,:,2), using
%
%   -kocc  also calculate occupation for each Wilson site
%          in global ground state; alternatively,
%   -kcum  also calculate cumulative occupation for sites up to site k1
%          in ground state space at iteration k
%   -v     verbose mode.
%
% Wb,Jul28,09 ; Wb,Dec14,10

% adapted from roncat.m
% Wb,Dec14,10 switched: sLR order => LRs order

  getopt('init',varargin);
     keep =getopt('-k');
     vflag=getopt('-v');
   % T   =getopt('T', 0);
     if     getopt('-kocc'), kocc=1; k1=1;
     elseif getopt('-kcum'), kocc=2;
     k1  =getopt('k1',4); else kocc=0; k1=-1; end
  nrg=getopt('get_last',[ getenv('LMA'), '/NRG/NRG' ]);

  load([nrg, '_info.mat']);

  L=length(ops.ff)+1; NN=NN(:); nops=numel(NN);
  if kocc<1
       nn=nan(L-1,nops,1);
  else nn=nan(L-1,nops,2); end

  RR=QSpace(L-1,2); X=[]; gotlast=0;

  nloc=NN;

% skip last iteration (creates large matrizes for 2-channel calculation
% check for MatLab memory management / fragmentation -- Wb,Apr27,09
  k2=L-1;

% RR(:,1) = density matrix of local ground state spaces
% RR(:,2) = density matrizes wrt. overall ground state
  for k=k2:-1:1
     if vflag, fprintf(1,'\r   get RK(%2g/%g) ... \r',k,N); end
     f=sprintf([nrg '_%02d.mat'],k-1);
     if k<k2, s=load(f,'AK'); else
        s=load(f,'AK','HK'); R=getrhoQS(s.HK);
     end
     RR(k,2)=R;
     R=contractQS(s.AK,[2 3],contractQS(s.AK,2,R,2,'conjA'),[3 2]); % LRs
  end

  for k=1:(L-1) % kloop
     if vflag, fprintf(1,'\r   get NN(%2g/%g) ... \r',k,N); end
     f=sprintf([nrg '_%02d.mat'],k-1);
     s=load(f,'AK','HK');

     if k==1
        if norm(diff(s.AK.Q{1},[],1)) % LRs order
        lflag=1; % got L-space for A0
        else lflag=0; end
     end

     if isempty(s.AK.Q), gotlast=1;
     s=load(f,'AT','HT'); s=struct('AK',s.AT,'HK',s.HT); end

   % nn(k,:,1) = dot occupation according to ground state at iteration k
     R=getrhoQS(s.HK); RR(k,1)=R;
     for j=1:nops
        if k>1 || lflag
             Q=contractQS(NN(j),2,s.AK,1); % keep LRs order
        else Q=contractQS(s.AK,3,nloc(j),2); end

        NN(j)=contractQS(s.AK,[1 3], Q, [1 3],'conjA');
        nn(k,j,1)=trace(R*NN(j));
     end

   if kocc==1
   % nn(k,:,2) contains occupation of each Wilson site
   % w.r.t. overall ground state; eg. AOC:
   % calculate displaced charge // tags: Wilson site occupation
     if lflag || k>1
     for j=1:nops
        Q=contractQS(s.AK,3,nloc(j),2); % keep LRs order
        Q=contractQS(s.AK,[1 3],Q,[1 3],'conjA');
        nn(k,j,2)=trace(RR(k,2)*Q);
     end
     end
   elseif kocc==2
   % nn(k,:,2) contains overall occupation of Wilson sites up to k1
   % similar to dot occupancy (ie. in ground state space of iteration k
     if k>k1
        for j=1:nops
           Q=contractQS(NN(j,2),2,s.AK,1); % keep LRs order
           NN(j,2)=contractQS(s.AK,[1 3],Q,[1 3],'conjA');
           nn(k,j,2)=trace(R*NN(j,2))-(k1-1+lflag)*0.5;
        end
     elseif k>1 || lflag
        for j=1:nops
           Q=contractQS(s.AK,3,nloc(j),2);
           if k>2 || k>1 && lflag
              Q=QSpace(Q)+contractQS(NN(j,2),2,s.AK,1); % keep LRs order
           end

           NN(j,2)=contractQS(s.AK,[1 3],Q,[1 3],'conjA');
         % nn(k,j,2)=trace(R*NN(j,2));
        end
     end
   end

     if gotlast, break; end
  end % of kloop

  if vflag, fprintf(1,'\r%30s\r',''); end

  if nargout>1, Iocc=add2struct('-',RR,NN,kocc,k1,k2); end

end

