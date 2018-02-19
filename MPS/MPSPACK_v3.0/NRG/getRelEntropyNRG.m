function [I,II,NRG,NRG0]=getRelEntropyNRG(NRG,NRG0,varargin)
% function getRelEntropyNRG(nrg,nrg0 [,opts])
% 
%    calculate entropy from NRG run <NRG> while comparing it
%    to another NRG run <NRG0> (eg. considered to be an NRG run
%    with passivated impurity).
% 
% Options
% 
%   '-p'    plot result
%   'TT',.. specify temperatures other than default (Wilson scales)
%   'z',..  apply z-shift to TT defined from Wilson scale
% 
% Wb,Aug12,11

% adapted from tst_entropy2.m

  if nargin<2 || xor(ischar(NRG),ischar(NRG0))
     helpthis, if nargin || nargout
     error('Wb:ERR','invalid usage'), end, return
  end

  getopt('init',varargin);
     TT    = getopt('TT',[]); if isempty(TT)
     z     = getopt('z',0); end
     pflag = getopt('-p');
  getopt('check_error');

  if ischar(NRG)
     i=sprintf('%s_info.mat',NRG);
     i0=sprintf('%s_info.mat',NRG0);
     if ~exist(i,'file') || ~exist(i0,'file'), { NRG, NRG0 }
        error('Wb:ERR','\n   ERR invalid NRG data (%s)',i); end
     Inrg=load(i);
  else
     Inrg=NRG;
  end

% ff=Inrg.param.ff;
  if isempty(TT)
   % TT=logspace(log10(ff(end)),2,256);
   % NB! stick with energy-shells => allows simple even/odd average below!
   % tags: EO_AVERAGE
     TT=Inrg.Lambda.^((-(Inrg.N-5)/2:0.5:4)+z);
  elseif z, wblog('WRN','z-shift (%g) will be ignored',z);
  end
  TT=reshape(TT,[],1);

  if ischar(NRG ), wblog(' * ','loading data %s',NRG ); end
  [II(1),NRG ]=specificHeatNRG(NRG, 'TT',TT);
  if ischar(NRG0), wblog(' * ','loading data %s',NRG0); end
  [II(2),NRG0]=specificHeatNRG(NRG0,'TT',TT);

  C=cat(2,II.C); Cimp=C(:,1)-C(:,2); cc=avgdata(Cimp); % EO_AVERAGE
  S=cat(2,II.S); Simp=S(:,1)-S(:,2); ss=avgdata(Simp); % EO_AVERAGE
  tt=sqrt(TT(1:end-1).*TT(2:end));

  I=add2struct('-',TT,Cimp,Simp,tt,cc,ss);

if ~pflag, return; end

  global param ; param=Inrg.param;

ah=smaxis(2,2,'tag',mfilename,'fpos',[890 470 710 650]);
header('%M :: %s (%s)',NRG.nrg,NRG0.nrg); nrg_header; addt2fig Wb

  h=header('fleft');
  set(h,'String',regexprep(get(h,'String'),'N_K=\d+',...
   sprintf('N_K{\\\\leq}%d (%d)',max(Inrg.NK(:,1:2)))));

setax(ah(1,1))
  bs=sprintf('FDM (ref)');

  semilogx(TT,C(:,1),'o-','Disp','FDM'); hold on
  semilogx(TT,C(:,2),'r','Disp',bs); hold on

  xtight; sms(2); label('T','total specific heat C');
  legdisp('Location','SouthWest');

setax(ah(1,2))

  semilogx(TT,S(:,1),'o-','Disp','FDM'); hold on
  semilogx(TT,S(:,2),'r','Disp',bs); hold on

  xtight; sms(2); label('T','total entropy S');
  legdisp('Location','SouthEast');

setax(ah(2,1))

  h=semilogx(TT,Cimp,'o-','Disp','FDM'); xtight; sms(2); blurl(h); hold on
  h=semilogx(tt,cc,'Disp','FDM (E/O-avg)'); % EO_AVERAGE
  label('T','relative C (impurity contribution)');

  if isfield(Inrg.param,'U')
       q=[Inrg.param.U, Inrg.param.Gamma, TKondo(Inrg.param)];
  else q=nan(1,3); end
  xmark(q);

  legdisp('Location','NorthWest');

setax(ah(2,2))

  tt=sqrt(TT(1:end-1).*TT(2:end));
  ss=avgdata(Simp);

if 0

  h=semilogx(TT,Simp/log(2),'o-','Disp','FDM'); hold on
  sms(2); blurl(h); xtight;

  h=semilogx(tt,ss/log(2),'Disp','FDM (E/O-avg)'); % EO_AVERAGE

  label('T','S_{imp}/ln(2)  (relative impurity contribution)');

  if 0
     ih=inset('NW','scale',[1.6 1],'dx',[0.1 0]);
        h=loglog(TT,exp(Simp),'o-'); hold on
        semilogx(tt,exp(ss)); blurl(h); sms(2); xtight;
        ytight(1.1); ymark(1:10,'k:'); grid on
        label('T','exp(S_{imp})');
     setax(ah(2,2))
  end

else
  h=loglog(TT,exp(Simp),'o-','Disp','FDM'); hold on
  sms(2); blurl(h); xtight;

  h=semilogx(tt,exp(ss),'Disp','FDM (E/O-avg)'); % EO_AVERAGE

  label('T','exp(S_{imp})  (relative impurity contribution)');

% for i=[1 2 4 8 16], ymark(log(i),'k:','istr',sprintf('ln(%g)',i)); end
  ytight(1.1); set(gca,'YTick',1:10); grid on

  if 0
     lo={'interpreter','latex','Color','k','FontW','normal'};
     h=ymark(sqrt(2),'k:','istr',{'$\sqrt{2}$',-0.07,lo{:}});
  end

end

  legdisp('Location','SouthEast');

end

