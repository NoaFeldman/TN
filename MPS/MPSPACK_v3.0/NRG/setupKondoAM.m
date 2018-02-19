
% setting up SU(3) symmetric Anderson/Kondo hybrid
% outsourced from rnrg_SUN.m // Wb,Apr28,13
% Wb,Sep26,11

  setdef('NC',2','Gamma',2,'JH',8); % LARGE_U
  % Gamma=2, JH= 8 => TK=~1E-08
  % Gamma=2, JH=10 => TK=~1E-10
  % Gamma=2, JH=12 => TK=~1E-14

  istr=sprintf('ph-symmetric Anderson/Kondo hybrid (%g-channel)',NC);

  setdef('Lambda',4','N',50,'Nkeep',2048,'Etrunc',6);

  if JH>1
     emax=max(10,0.30*NC*JH*(NC*JH+2));
     if isset('Etrunc') && JH>Etrunc && ~isset('ETRUNC')
        if NC==1
         % start truncation with iteration NRG_05
         % chopping off the high-energy band due to large JH
         % n=round(2*log(4*JH/Etrunc)/log(Lambda)); % n=4 for Lambda=4
           if Lambda<2, n=6; else n=5; end
           setdef('NEX',n);
           ETRUNC=[ repmat(max(50*JH,Etrunc),1,NEX), JH ];
        else
           ETRUNC=repmat(max(100*JH,2*Etrunc),1,3);
        end
     end
  elseif JH<0
     wblog('WRN','got negative JH (%g)',JH);
  end

  global param g_TKondo

  co={}; if isset('z'), co(end+1:end+2)={'z',z}; end
  if exist('ALambda','var') && ~isempty(ALambda) 
     switch ALambda
        case { 1,'-AL'}, co{end+1}='-AL';
        case {-1,'-w' }, co{end+1}='-w';
        otherwise error('Wb:ERR','\n   ERR invalid ALambda');
     end
  end
  ff=getNRGcoupling(Gamma,Lambda,N,co{:});

  if isset('SYM'), sym=SYM; else
     if NC>1
        if isset('B')
             sym='Acharge,Aspin,SUNchannel';
        else sym='Acharge,SU2spin,SUNchannel'; end
     else
        if isset('B')
             sym='Aspin,SU2charge'; %% ERR
        else sym='SU2spin,SU2charge'; end
     end
  end

  param=add2struct('-','istr?',sym,NC,Gamma,JH,'B?',Lambda,N,Nkeep,...
    'z?','NKEEP?','Etrunc?','ETRUNC?',ff,'ALambda?');

% expected TKondo (see tex notes)
  if NC==1
     param.iTKondo='NB! equivalent SIAM parameters for TKondo';
     param.U=(3/2)*JH;
     param.epsd=-param.U/2; clear global g_TKondo
  elseif isempty(g_TKondo)
   % g_TKondo=sqrt(JH*Gamma/2)*exp(-JH/(2*Gamma) * (1-1/(2*NC+1)));
     g_TKondo=sqrt(JH*Gamma/2)*exp(-0.5*(NC+1)*JH/Gamma); % 4=2S+1 = 2(S+1/2)
     % factor (1/2) in exponent is empirical
  end

  disp(param);

  [FF,Z,S3,IS]=getLocalSpace('FermionS',sym,'NC',NC,'-v');

  if all(isAbelian(FF)), F3=FF; p3=[];
     FF=squeeze(FF);
     S3=squeeze(S3);
  else
     if NC==1, ZFLAG=2; end
     p3=[1 3 2];
  end

  S2=QSpace;
  for i=1:numel(S3)
     if numel(S3(i).Q)>2, ic=[1 3]; else ic=1; end
     S2=S2+QSpace(contractQS(S3(i),ic,S3(i),ic,'conjA'));
  end

  H0=-JH*S2; % i.e. ferromagnetic Hund's coupling for JH>0

% operators for self-energy: FN = comm(H0,FC)
  FN=FF;
  
  for i=1:numel(FN)
     FN(i)=QSpace(contractQS(H0,2,FN(i),1))-contractQS(FN(i),2,H0,1,p3);
     FN(i).info.otype=FF(i).info.otype;
     FN(i)=(1/JH)*FN(i); %  eliminate JH dependence (since H0 \ propto JH)
  end

if isset('init_only')
  if numel(FF(1).Q)==3 && dim(FF(1),3,'-f')
     for i=1:numel(FF)
        FF(i).Q(3)=[];
        FF(i).info={};
     end
  end

  A1=QSpace(getIdentityQS(FF,FF));
  A2=QSpace(getIdentityQS(A1,3,FF)); %  260 x 10 x   9086 =>   4096 x 64 x    262.144
% A3=QSpace(getIdentityQS(A2,3,FF)); % 9086 x 10 x 366744 => 262144 x 64 x 16.777.216
% A2 takes about <2 seconds
% A3 takes about 12 seconds
end

% generate auxiliary operator for empty (=vacuum) state to generate A0
% Q=H0; Q.Q=repmat({zeros(1,size(A1.Q{1},2))},1,2);
% Q.data={1}; Q.info.cgs=repmat({1},1,size(H0.info.cgs,2));
% A0=QSpace(getIdentityQS(Q,FF));

% A0=permute(addsingleton(IS.E,'-p'),[1 3 2]);
  A0=getAtensorLoc(IS.E); % Wb,Dec28,14

% operator setup of fdmNRG

  setdef('seq',0);
  if exist('NKEEP','var'), n=max([NKEEP Nkeep]); else n=Nkeep; end

  setdef('wops','FNS'); clear ac_cmd cflags

  switch wops
    case 'Sz'
       op1=[]; op2=S3(end); zflags=[0];
    case 'F1'
       op1=[]; op2=FF(1); zflags=[0];
    case 'FNS'
       if n<1800 || seq || NC<3
        % op1=[]; op2=[FF,S3]; zflags=[1 0]; % kinit
          op1=[FF(1),FN(1),S3(end)];
          op2=[FF(1),FF(1),S3(end)]; zflags=[1 1 0];

        % also calculate improved spectral function
          ac_cmd='[ac,Ic]=getGC_JJH(ox,ax(:,1),-JH*ax(:,2),''-q'');';
       else
          op1=[]; op2=FF; zflags=1;
       end
    otherwise error('Wb:ERR','\n   ERR invalid switch (wops)');
  end

