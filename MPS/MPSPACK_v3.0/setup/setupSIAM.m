% -------------------------------------------------------------------- %
% setup general parameters
% -------------------------------------------------------------------- %

  global param

  setdef('Lambda',2); setdef('N',ceil(50*log(2)/log(Lambda)));
  if ~exist('U','var')
       U=0.12; epsd=-U/3; Gamma=0.01; B=0;
  else initdef('B',0); end

  param=struct( ...
     'U', U, 'epsd', epsd, 'Gamma', Gamma, 'B', B, ...
     'N', N, 'Lambda', Lambda  ...
  );

  if exist('TK0','var')
  TK=TK0; else
  TK=TKondo; end

  if isstr(B), B=eval(B); param.B=B; end

% pg_r = pseudo-gag coefficient r
  co={'-x'};
  if exist('z',   'var'), co(end+1:end+2)={'z',z   }; param.z=z;    end
  if exist('pg_r','var'), co(end+1:end+2)={'r',pg_r}; param.r=pg_r; end

  if exist('ALambda','var')
     if ischar(ALambda) || isempty(ALambda) || ALambda
        co(end+1)={'-AL'}; param.ALambda=1;
     else
        co(end+1)={'-w'}; % param.ALambda=0;
        if isfield(param,'ALambda') param=rmfield(param,'ALambda'); end
     end
  end

  if ~exist('FFLANCZOS','var')
     [ff,N,Ic]=getNRGcoupling(Gamma(1),Lambda,N,co{:}); param.ff=ff;
  else
     if ~isvector(FFLANCZOS), error('Wb:ERR',...
        'setupSIAM - FFLANCZOS must be a vector!'); end
     m=1+length(FFLANCZOS);

     if N>m,     s=sprintf('N=%d->%d',N,m); N=m;
     elseif N<m, s=sprintf('N=%d/%d',N,m);
     else        s=sprintf('N=%d',N); end

     wblog('NB!','Using couplings FFLANCZOS (%s)',s);

     ff=FFLANCZOS(1:N-1);
     param.N=N;
  end

% allow string expressions for B and BX containing `TK'
  if exist('B','var') && ischar(B)
  eval(['B=' B ';']); param.B=B; end
  if exist('BX','var') && ischar(BX)
  eval(['BX=' BX ';']); param.BX=BX; end

  if exist('BX','var') && isscalar(BX)
     if exist('QDIM') && ~isequal(QDIM,1)
        wblog('WRN','Since BX=%.3g, set QDIM=%d -> 1',BX,QDIM);
     end
     QDIM=1;
  end

  wblog('<i>', 'SIAM parameters (TK=%.4e)\N', TK);
  disp(param);

% -------------------------------------------------------------------- %
% operator setup
% -------------------------------------------------------------------- %

  if ~exist('QDIM'), QDIM=2; end
  if isequal(QDIM,2)

   % Hamiltonian H0 and its basis A0
   %    [ charge, 2*Sz ] basis // up down
     QS = [-1  0              %    0  0    empty
            0 -1              %    0  1    spin down
            0 +1              %    1  0    spin up
           +1  0 ];           %    1  1    doubly occupied

            %    spin down,spin up => H=-BSz (spin up has lower energy for B>0)
     data = [ 0, epsd+B/2, epsd-B/2, 2*epsd + U ] + U/2;

     H0 = QSpace( ...
            { QS, QS }, ...
            mat2cell( data, 1, ones(length(data),1) )' ...
          );

     A0=QSpace(zeros(1,size(QS,2)), QS, 'identity');

   % coupling operators
     clear FC op1 op2

     FC(1)=QSpace(...  % spin up
             [-1  0; 0  1],  1, ...
             [ 0 -1; 1  0],  1  ...
           );

     FC(2)=QSpace(...  % spin down
             [-1  0; 0 -1],  1, ...
             [ 0  1; 1  0], -1  ...
           );

     FN(1)=QSpace([ 0 -1; 1  0],  1);
     FN(2)=QSpace([ 0  1; 1  0], -1);

     Z = QSpace(...
             { QS, QS }, ...
             mat2cell( (-1).^[0 1 1 2], 1, ones(4,1))' ...
          );

     SZ = QSpace(...
        [ 0  1; 0  1], +1, ...
        [ 0 -1; 0 -1], -1  ...
     );

   % local operators

     if ~exist('gg','var') || isempty(gg)
        gg=[]; FL=QSpace;
     else
        if numel(gg)==1, gg=repmat(gg,size(ff));
        elseif numel(gg)<numel(ff), error('Wb:ERR',['\n   ' ...
          'ERR invalid dimension of gg (local couplings)']); 
        end
      % eg. use FL=N0 below
        if ~exist('FL','var') || min(size(gg))~=numel(FL)
           if ~exist('FL','var'), n=0; else n=numel(FL); end
           error('Wb:ERR',['\n   ERR size mismatch for ' ... 
           'local operators (%g/%g)'],min(size(gg),n));
        end
        wblog(' * ','using local operators');
     end

   % operators for dmNRG
     Z0=Z;

     if exist('isBulla')==1 && isBulla
        C = [FC FC];
        F = [FC FC];

        NN(1) = QSpace(QS, diag([0 1 0 1]), 'operator'); % = number n2 (!)
        NN(2) = QSpace(QS, diag([0 0 1 1]), 'operator'); % = number n1 (!)

        for i=1:2
        F(i+2)=skipzeros(QSpace(contractQS(F(i+2),2,NN(i),1))); end
     else
        F=[]; C=FC;

        if isset('B') || isset('Bflag') % Wb,Sep28,11
           op1=[]; op2=[FC, SZ]; zflags=ones(size(op2)); zflags(end)=0;
        else
         % include magn. suszeptibility to determine TKondo
           op1 = [ FC(1) FN(1) SZ ];
           op2 = [ FC([1 1]) SZ ]; zflags=[1 1 0];
         % used in dma_plot.m
           ac_cmd='[ac,Ic]=getGC_bulla(ox, ax(:,1),  ax(:,2));';
        end

     end

     clear NN i % QS

% -------------------------------------------------------------------- %
  elseif isequal(QDIM,1)
% -------------------------------------------------------------------- %

   % setup initial Hamiltonian: QSpace H0 and first A matrix
   %   { charge } basis; up down
     QS = [ 0          %  0  0
            1          %  1  0
            1          %  0  1
            2 ];       %  1  1

     %    NB! using H=-BSz (spin up has lower energy for B>0)
     H0 = diag([ 0, epsd-B/2, epsd+B/2, 2*epsd + U ] + U/2);
     if exist('BX','var') && isscalar(BX)
        H0(2,3)=BX/2;
        H0(3,2)=BX/2;
     end

     H0=QSpace(QS,H0,'operator');
     A0=QSpace(0,QS,'identity');
     Z =QSpace(QS,diag((-1).^QS),'operator');

   % setup local operators: row vector<QSpace> FL(m)

     FL=QSpace; gg=[];

   % setup coupling operators

     FC(1)=QSpace(QS, [  % spin up
            0  0  1  0
            0  0  0  1
            0  0  0  0
            0  0  0  0 ],'operator');

     FC(2)=QSpace(QS, [  % spin down
            0  1  0  0
            0  0  0  0
            0  0  0 -1
            0  0  0  0 ],'operator');

   % operators for dmNRG
     Z0=Z; F=[]; C=FC;

     clear NN i % QS

% -------------------------------------------------------------------- %
  else error('Wb:SIAM','Invalid QDIM'); end
% -------------------------------------------------------------------- %

  E0=QSpace(contractQS(A0,[1 2],A0,[1 2]));

  N0=FC;
  for i=1:numel(FC)
     N0(i)=contractQS(FC(i),1,FC(i),1);
  end

