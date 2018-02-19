function varargout=getEDdata(varargin)
% function [[EK,]ED,Inrg]=getEDdata([NRG ,opts])
%
%    get energy data of discarded state space from specified NRG run
%    (default: 'NRG/NRG'). All energies are converted to physical
%    energies w.r.t. the ground state energy of the last iteraion,
%    i.e. all NRG scaling and shifting are undone.
%
%    The returned cell array ED{k} contains one matrix [Ed,deg]
%    for every Wilson iteration k.
%
% To be used for calculating specific heat, entropy, and similar.
%
% Options
%
%    -K   also return combined data for each iteration (i.e. K+D) as first argument
%    -Q   include symmetry labels
%    -m   merge ED data into a single 3-colmun matrix [k,Ed,deg]
%         of Wilson iteration k, discarded energies Ed and
%         their respective degeneracies
%    -d   get also local dimensions along Wilson chain from AK
%         (returned as Inrg.d)
%    -q   quiet mode.
%
% Wb,Jul17,11

% adapted from specificHeatNRG.m

  getopt('INIT',varargin);
     merge=getopt('-m');
     Kflag=getopt('-K');
     dflag=getopt('-d');
     qflag=getopt('-q');
     Qflag=getopt('-Q');
  NRG=getopt('get_last','NRG/NRG');

  m=[NRG '_info.mat']; mf=mfilename;
    if ~exist(m,'file'), NRG
        error('Wb:ERR','\n   ERR invalid NRG data'); end
    if ~qflag, fprintf(1,'\r   %s: loading %s ...  \r',mf,repHome(m)); end
  Inrg=load(m);
    if ~qflag, fprintf(1,'\r%80s\r',''); end

% ES=0.5*(Lambda+1).*Lambda.^(-0.5*(0:length(Inrg.E0)-1)); % energy scale
  ES=Inrg.EScale;
  E0=Inrg.E0.*ES;

% dE=cumsum(E0); Eref=min(dE); dE=dE-Eref;
% NB! this looses energy resolution towards the end of the system
% NB! with Eref=sum(E0), above cumsum() is equivalent to performing
% (-1)*cumsum() wrt. to the end of the system! => numerically stable.
% Wb,Jul19,13

  dE=-fliplr(cumsum(fliplr(E0)));
  dE=[ dE(2:end), 0 ];

  if min(dE)<0, wblog('WRN','got min(dE)=%.3g !??',min(dE)); end
  if abs(Inrg.phE0-sum(E0))/Inrg.phE0>1E-8, error('Wb:ERR',...
    'ground state energy discrepancy (%.3g)',sum(E0)-Inrg.phE0);
  end

  if dflag, vars={'AK','HT'}; else vars={'HT'}; end
  if Kflag, vars{end+1}='HK'; end

  k=0; HT={};
  while 1
     m=sprintf('%s_%02d.mat',NRG,k);
     if ~exist(m,'file'), break; end
     if ~qflag, fprintf(1,'\r   %s: loading %s ...  \r',mf,repHome(m)); end

     k=k+1;
     load(m,vars{:});

     gotHT=(~isempty(HT) && ~isempty(HT.data));

     if k==1 % Wb,Feb26,14
        if Kflag
           if isdiag(QSpace(HK),'-d')<2
              [~,q]=eigQS(HK); AK=q.AK; HK=q.EK;
           end
        end
        if gotHT, error('Wb:ERR',...
          '\n   ERR got truncation at 0th NRG iteration !??');
        end
     elseif k==2
        if Kflag
           if isdiag(QSpace(HK),'-d')<2, error('Wb:ERR',...
             '\n   ERR got non-diagonal HK at 1st NRG iteration !??');
           end
        end
     end

     if gotHT
      % physical energy sets
        ED(k).E=ES(k)*cat(2,HT.data{:})' + dE(k);
        ED(k).deg=getzdim(QSpace(HT),2,'-p','-x'); % degeneracies

        if Qflag
           nd=numel(HT.data); qq=cell(1,nd); Q=HT.Q{1};
           for i=1:nd
              qq{i}=repmat(Q(i,:),length(HT.data{i}),1);
           end
           ED(k).Q=cat(1,qq{:});
        end
     end

   % consider only iterations where truncation occured
   % eg. allow impurity to have altered state space
     if dflag && ~isempty(HT) && ~isempty(HT.data) ...
              && ~isempty(AK) && ~isempty(AK.data)
        q=getDimQS(AK); d(k)=q(end);
     end

     if Kflag
        if isempty(HT.data) && ~iscell(HT.data), HT.data={}; end
        if isempty(HK.data) && ~iscell(HK.data), HK.data={}; end

        EK(k).E=ES(k)*[cat(2,HK.data{:}), cat(2,HT.data{:})]'; % + dE(k);
        EK(k).deg=[ % degeneracies
            getzdim(QSpace(HK),2,'-p','-x')
            getzdim(QSpace(HT),2,'-p','-x') ];
        if Qflag
           nd=numel(HK.data); qq=cell(1,nd); if nd, Q=HK.Q{1}; end
           for i=1:nd
              qq{i}=repmat(Q(i,:),length(HK.data{i}),1);
           end
           qq=cat(1,qq{:}); if gotHT
           EK(k).Q=[qq; ED(k).Q ]; else EK(k).Q=qq; end
        end

      % sort data w.r.t. energy
        [EK(k).E,is]=sort(EK(k).E);
        EK(k).deg=EK(k).deg(is); if Qflag
        EK(k).Q  =EK(k).Q(is,:); end
     end
  end
  if ~qflag, fprintf(1,'\r%70s\r',''); end

  if dflag, 
     Inrg.d=d;
     d=unique(d(find(d))); if numel(d)~=1, d
       error('Wb:ERR','\n   ERR failed to determine local d'); end
     Inrg.dloc=d;
  end

  if k<6, wblog('WRN','Wilson chain appears to short'); end
  N=k;

  if merge
     for k=1:N
        if isempty(ED(k).E), continue; end
      % combined ED data into single matrix to ED = [ k, E, deg [,Q] ]
        q={ repmat(k-1,numel(ED(k).E),1), ED(k).E, ED(k).deg };
        if Qflag, q{end+1}=ED(k).Q; end
        ED(k).dd=cat(2,q{:});
     end
     ED=cat(1,ED.dd);
  end

  if ~nargout, varargout={ED}; else
     if Kflag
          varargout={EK,ED,Inrg};
     else varargout={ED,Inrg}; end
     varargout=varargout(1:nargout);
  end

end

