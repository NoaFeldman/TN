function [dw,Iout,IR]=getDiscWeightNRG(varargin)
% function [dw,I,IR]=getDiscWeightNRG(NRGtag [,opts])
%
%    get estimate for discarded weight (dw) from reduced density
%    matrices, considering backward update from overall ground
%    state space (default NRGtag: './NRG/NRG').
%
% Options
%
%   'n0',...   how many times to back-propagate rho until min(rho) is obtained
%   'k0',...   consider first NRG_<it> iterations with it<k0 exact
%   'beta',..  local 1/T to initialized rho to ground state space (100)
%   'chi',..   lower range of eig(rho) to use (default: 0.05 = last 5%)
%
%   '-E'       use energy eigenbasis (E-basis s, rather than
%              eigenbasis of rho, R-basis r)
%   'Etr',..   truncation energy (as used in histogram for energy)
%              (default: Inrg.Itr.Etrunc(1)) 
%
%   '-q'       quiet (do not show iteration k while calculating rho)
%   '-fig'     show summarizing figure
%
% The returned (estimated) discarded weight is given by dw = mean(I.e1).
% where the returned info structure I contains (amongst others)
% the following fields:
%
%   .chi   lowest range of eigenvalues used (0.05 = 5%)
%   .e1    discarded weight based on the lowest weight states
%          in the reduced density matrices (RDM) R(n-n0:n);
%          (lowest chi percent in terms of number of states)
%   .e2    discarded weight based on the top energy states
%          in the reduced density matrices (RDM) R(n-n0:n);
%          however, since rho_r and E_r are strongly related,
%          e2 is typically only slightly larger than e1 (R-basis).
%   .dw    estimated discarded weight for each Wilson shell
%   .se    Shannon entropy
%   .beta  1/T (in rescaled units) used when initializing the density
%          matrix in the kept energy states (projection into
%          the ground state space; default: beta=100).
%
% Wb,Jan24,11

% adapted from: getrhoNRG_resall.m ; Wb,Aug17,12

  getopt('init',varargin);

     n0    =getopt('n0',-1);
     k0    =getopt('k0', 0);
     beta  =getopt('beta',200);

   % smooth=getopt('smooth',4); % Wb,May12,14
   % NB! estimate for disc. weight will always be discontinuous in Lambda
   % since energy truncation is sharp; smoothing contribution to e1 or e2
   % below cannot change this fact // Wb,May12,14
     smooth=getopt('smooth',0);

   % take top 10% in energy to evaluate res_disc
     if smooth<=0
        chi =getopt('chi',0.05);
     else
      % NB! using smooth, takes sum() instead of mean() of
      % low-energy spectrum; together with smoothing, this
      % results in an effectively higher chi
      % see SMOOTH below // Wb,May12,14
        chi =getopt('chi',0.01);
     end

     vflag =getopt('-v');    % write final log entry for discarded weight
     qflag =getopt('-q');     % do not show progress (current iteration)
     Eflag =getopt('-E');     % Wb,Nov19,15
     Etr   =getopt('Etr',[]);  % Wb,Nov20,15
     showfig=getopt('-fig'); % plot figure

     yl    =getopt('yl',[]); % ylim range for epsilon^D_n panel ah{2,1}(end)

     aflag =getopt('-a');
     xl31  =getopt('xl31',[]);
     yl31  =getopt('yl31',[]);

     nochi =getopt('~chi');  % do not plot chi (affects plot only)
     nx    =getopt('nx',[]);

     fflag =getopt('-ff');

  varargin=getopt('get_remaining'); narg=length(varargin);

  if narg
     if narg>1 || ~ischar(varargin{1}), varargin
     error('Wb:ERR','invalid NRGdata tag'); end
     nrg=varargin{1};
  else nrg='./NRG/NRG'; end

  if numel(chi)~=1 || (chi<=0 || chi>=1)
     error('Wb:ERR','\n   ERR invalid chi (%g)',chi);
  end

  if smooth>0 && smooth<2, error('Wb:ERR',... % tags: SMOOTH
    '\n   ERR invalid smoothing parameter for disc. weight (%g)',smooth);
  end

  if nrg(1)~='/' && isempty(findstr(pwd,'Data')), cto lma, end
% ff=dir2([nrg '_\d+.mat']);
  ff=getnrgfiles(nrg);
% ensures that file names are sorted => data is in Wilson shell order
  N=numel(ff); 

  if ~N, error('Wb:ERR',...
  '\n   ERR no NRG files ''%s'' found',nrg); end

% NB! *_info also contains ff(!), and overwrites global param if loaded!
  Inrg=load([nrg '_info']);
  Lambda=Inrg.Lambda; Nkeep=Inrg.Nkeep;

  if isempty(Etr)
     Etr=Inrg.Itr.Etrunc(1); % Etrunc(2) is relEtrunc if any
  end
  if numel(Etr)~=1 || Etr<=0, Etr
     error('Wb:ERR','\n   ERR invalid usage (Etr)');
  end

  mat=ff(N).name; load(mat);
     if ~isempty(QSpace(HK)), error('Wb:ERR',...
     'HK must be empty at last iteration (%s) !??',mat); end

  AK=AT; AT=QSpace;
  HK=HT; HT=QSpace;

  if fflag, getbase FC
      if isempty(FC), wblog('ERR',...
        'failed to find FC operators (got -ff flag)'); fflag=0;
      else
         for i=1:numel(FC)
          % NB! trace(FF') is a scalar! => plain weight factors x
            x(i)=normQS(FC(i))^2;
         end
      end
   end
   e3=[];

 % -------------------------------------------------------------------- %
 % build reduced density matrix space starting from ground state space
  if aflag
     if isequal(n0,-1), n0=[]; end
     n1=1; n2=N;
  else
     if isequal(n0,-1), % n0=4:8;
      % adjust n0 with Lambda to have Lambda^(n0/2) = const
      % with lower bound n0>=[3 4] // Wb,Nov20,15
        q=max(2, 4/log(Lambda));
        % q=9.87 for Lambda=1.5
        % q=5.77 for Lambda=2
        % q=2.88 for Lambda=4
        % q=1.90 for Lambda=8 -> q=2

      % ensure that d^n0 > Nkept
        d=getDimQS(Inrg.ops.A0); d=d(end);
        k=max(Inrg.NK); if numel(k==4) k=k(2); else k=k(1); end
        q=max(q, log(k)/log(d));

        n0=floor(q)+[1 2]; % average at least over even/odd
     end

     if ~isempty(n0), n1=min(n0); n2=max(n0);
     else, n1=4; n2=4; % n0 set empty by input
     end
  end

% estimate discarded weight for each shell along the Wilson chain
  for k=N:-1:2

     [R,I0]=getrhoQS(HK,beta); % also calls skipzeros(R)
     IR(k)=struct('k1',[],'k2',k,'R',R,'R0',R,'I0',I0); % Rho(k1:k2)

   % t=trace(R); if abs(t-1)>1E-12, t,t-1, keyboard, end

   % error estimator trace( rho V_KD (V_KD)' )
     if fflag && ~isempty(HT) && ~isempty(HT.Q), for i=1:numel(FC)
      % take correction to ground state due to 2nd order PT
      % [see notes: "fidelity susceptibility"] // Wb,Sep30,11
      % take sqrt() => allows to create mirror object A to be
      % contracted onto itself (line [a])
        R=getrhoQS(HK,beta,'-sqrt');
        q=contractQS(contractQS(AK,2,R,1),2,FC(i),1); % LRs(o)
        q=contractQS(q,'13',AT,'13*'); %  R(o),R'
        l=numel(q.Q);

        a=HT; for j=1:numel(HT.data), a.data{j}=1./a.data{j}; end
        q=contractQS(q,l,diag(QSpace(a)),1);
      % ic=1:length(q.Q);
      % q=getscalar(QSpace(contractQS(q,ic,q,ic))); % line [a]
        q=normQS(q)^2; % line [a] (same as above, in short)

        e3(k,i)=( x(i)*q*Inrg.ops.ff(k)^2 )^2;
        % 2nd order perturbative result to change in ground state wave function
        % => corresponds to singular values! => take (...)^2
     end, end

     if k>1
     load(ff(k-1).name,'HK','HT'); end


   % kk: which IR(kk) need to be backpropagated
   % (always include ground state space from last iteration)
     kk=unique([ (k-1)+(1:n2), N ]); kk(find(kk>N))=[];

     if ~qflag, fprintf(1,...
       '\r   %s %4g/%g (%g)... \r',mfilename,k,N,numel(kk)); end

     for k2=kk
        R=IR(k2).R;

        q=contractQS(AK,'2*',IR(k2).R,2); % LRs order (!)
        IR(k2).R=QSpace(contractQS(AK,'23',q,'32'));  % LRs order (!)
        IR(k2).k1=k-1;

        if k2-k+1<n1
         % haven't reached yet n0 back propagations for Rho_reduced
           continue
        end

        s.I0=IR(k2).I0;
        s.k1=IR(k2).k1;
        s.k2=IR(k2).k2; if k2==N, RR(k)=IR(k2); end

        [rr,Ie]=eigQS(IR(k2).R); q=Ie.EK; n=numel(q.data);
        for j=1:n, d=q.data{j}; d(find(d<1E-16))=nan;
            q.data{j}=fliplr(-log(d));
        end
        s.ES=q;

      % s.se=SEntropy(s.rr,'-rho'); this ignores CGC degeneracy!
        [s.se,s.rr,s.Is]=SEntropy(QSpace(IR(k2).R),'-xm');
      % NB! rr = eig(rho) expanded in multiplet dimension => sum(rr)=1.
      % NB! rr is already sorted in descending order

      % E1_chiN :: last chi percent of *ALL STATES*
      % (multiplets are already expanded)
        l=ceil((1-chi)*length(s.rr));
        if smooth>1
         % avoid sharp transition in weights which would lead to
         % discontinuities of dw as function of Lambda, etc.
           s.e1=(1./((s.rr/s.rr(l)).^smooth+1))*s.rr';
         % this is the Fermi function on a logscale using beta:=smooth(=4)
         % [with the variable beta already used in a different context above]
         % (NB! smooth>1 required to suppress contribution of rr>>rr(l) data)
         % In contrast to mean() below, here e1 is the SUM of all discarded
         % weights together with a smooth transition weighting function.
         % Wb,May12,14 // tags: SMOOTH
        else
           ix=find(s.rr<=1.001*s.rr(l)); % do not split degenerate spaces
           s.e1=mean(s.rr(ix));
        end

        if s.e1<0
           if abs(s.e1)>1E-12, wblog('ERR','got s.e1=%g',s.e1); end
         % s.e1=0; % keyboard
        end

      % match top chi percent wrt. to energy E=<r|H|r> % E2_chiE
      % as compared to just the smallest chi% eigenvalues of R.
      % => get <r|H|r> vs. eig(rho) data
        if isempty(HT.data)
         % no truncation => no discarded weight!
           s.erd=[]; s.e1=[]; s.e2=[];
        else
           R=struct(IR(k2).R);
           [i1,i2,I]=matchIndex(R.Q{1},HK.Q{1});
           dz=getzdim(QSpace(R),2,'-p');

           for i=1:numel(i1)
              j1=i1(i); ri=R.data{j1};
              j2=i2(i); ei=HK.data{j2};

            % do not include multiplet dimensions dz in eig(rho) weights
            % => need dz in averaging below! (mean) // Wb,May10,14
              di=repmat(dz(j1),size(ri,1),1);

              if Eflag
               % use energy eigenbasis (s) // Wb,Nov19,15
                 ri=diag(ri); ei=ei';
              else
               % using eigenbasis of R (r)
                 [u,ri]=eig(ri+ri'); ri=0.5*diag(ri);
                 ei=diag(u'*diag(ei)*u);
              end
              R.data{j1}=[ei,ri,di];
           end
           s.erd=cat(1,R.data{:});

         % q=[ min(s.erd(:,1)), max(s.erd(:,1)) ];
           ex=(1-chi)*max(s.erd(:,1));

           if smooth>1
              s.e2=(1./((s.erd(:,1)/ex).^smooth+1))'*(s.erd(:,2).*s.erd(:,3));
            % NB! in contrast to e1 above (which uses smoothing in r),
            % e2 refers to the last chi percent in E and hence also
            % uses smoothing in E! // Wb,May12,14 // tags: SMOOTH
           else
            % NB! do not split (multiplet) degeneracies
              ix=find(s.erd(:,1)>ex);

            % if numel(ix)<2, wblog('WRN',['got only %g/%g points ' ...
            %   '(chi=%.3g, k=%g)'],numel(ix),size(s.erd,1),chi,k);
            % end

            % s.e2=mean(s.erd(ix,2)); // Wb,May10,14
              s.e2=(s.erd(ix,2)'*s.erd(ix,3))/sum(s.erd(ix,3));;
           end

        end

      % diagonal elements => energy eigenbasis
        r=struct(diag(IR(k2).R));
        [s.rd,i]=sort(cat(1,r.data{:})); % rd = diagonal diag(rho)
        s.ed=min(s.rd);
      % l=length(s.rd); s.e2=mean(s.rd(1:ceil(l*chi)));
      % s.e2=0; % see below

      % II is N x (numel(n0)+1) structur array, where the last column
      % contains discarded weight derived from overall ground state.
        if k2-k<n2, II(k-1,k2-k+1)=s; end
        if k2==N, II(k-1,n2+1)=s; end % store data related to last iteration

        if abs(trace(IR(k2).R)-1)>1E-8, IR(k2).R % safeguard
        error('Wb:ERR','invalid density matrix R'); end
     end

     if k>1
        if fflag, v={'AT'}; else v={}; end
        load(ff(k-1).name,'AK',v{:},'HT');
     end

  end
  if ~qflag, fprintf(1,'\r%50s\r',''); end

  Iout=add2struct('-',nrg,Nkeep,Etr,beta,smooth,e3,n0,chi,fflag,aflag);
  if aflag && nargout<3, Iout.IA=II; end

  II=II(:,[n0 end]);

  Iout.param=Inrg.param;

% NB! e1 may be negative due to numerical noise when computing eig(rho)!
% Wb,Jan09,14
  se=getdatafield(II,'se','-q'); se=seteps2nan(se); Iout.se=se;
  e1=getdatafield(II,'e1','-q'); [e1,gotneg_e1]=seteps2nan(e1);
  e2=getdatafield(II,'e2','-q'); e2=seteps2nan(e2);

% -------------------------------------------------------------------- %
% NB! e2 is preferential in particular also when using Eflag
% due to the large spread of lowest-weight states!
% NB! Eflag leads to markedly smooter data when plotting
% discarded weight vs. Wilson shell!
% Wb,Nov19,15
% -------------------------------------------------------------------- %

% NB! both e1 and e2 are in eig(rho) basis, yet are calculated slightly
% differently (last chi percent in full set / wrt. energy, respectively.
  e=e1(k0+1:end,:); e(find(isnan(e)))=0; Iout.e1=max(e, [],1); % E1_chiN
  e=e2(k0+1:end,:); e(find(isnan(e)))=0; Iout.e2=max(e, [],1); % E2_chiE

% dw_out is mean across n0 (max otherwise!)
  dw=mean(Iout.e2); Iout.dw=max(e2,[],2); % E2_chiE (as defined in paper)
% dw=mean(Iout.e1); Iout.dw=max(e1,[],2); % E1_chiN

  if gotneg_e1~=0 && (dw<=0 || ~qflag || vflag)
     wblog('WRN','got e1<0 !! (%.3g%%)',100*gotneg_e1);
  end

  if vflag
  wblog('<i>','max. discarded weight: %.3g',dw); end

% -------------------------------------------------------------------- %
  if ~showfig, return; end
% -------------------------------------------------------------------- %

  EE=Inrg.EE; nE=min(512,size(EE,1)); % number of energies to show

% recalculate / check n0 in middle of chain
  k=round(size(II,1)/2);
  k1=cat(2,II(k,:).k1); if norm(diff(k1))>1E-12
     error('Wb:ERR','inconsistent k1''s !??'); end
  k2=cat(2,II(k,:).k2);
  n0_=k2(1)-k1(1);

  t=mfilename;

  if Eflag, bs='_s'; else bs='_r'; end

ah=smaxis(8,1,'tag',t,'dx',0.10,'dy',0.007,'DY',0.02); 

  s=repmat({''},1,4); getbase param

  if isfield(param,'istr'),s{1}=param.istr; end
  if isfield(param,'sym'),s{2}=param.sym;
  elseif isfield(param,'wnrg'),s{2}=param.wnrg; end

  if k0>0, s{3}=sprintf('k_0=%g, ',k0+1); end
  if Eflag
       s{4}='{\color{red}using energy eigenbasis}';
  else s{4}='using eigenbasis of \rho';
  end

  evalin('caller','nrg_header'); header(['%M :: %s (%s)' 10 ...
    '\\chi=%.3g%% (%sn_0=%g, m=%g, \\beta=%g, smooth=%g, %s)'],...
     s{1:2},100*chi,s{3},max(n0),numel(n0),beta,smooth,s{4});
  set(gcf,'DefaultAxesLineW',1.5); addt2fig Wb

  ah=mat2cell(ah,[2,3,3],1); dy=0.012;

  mvaxis(ah{1}(1),[0  dy+.003]);
  mvaxis(ah{1}(2),[0  dy+.006]);
  mvaxis(ah{2}(1),[0  dy+.006]);
  mvaxis(ah{2}(2),[0  dy+.009]);
  mvaxis(ah{2}(3),[0  dy+.012]);

  mvaxis(ah{3},[0  -.006]);

  p=get(ah{3,1}(3),'Pos'); p(3)=0.6*p(3); p(4)=2.2*p(4);
  set(ah{3,1}(1),'Position',p);

% ah{3,1}(2): top panel to (3,1)
  set(ah{3,1}(2),'Position',[p(1), p(2)+p(4)+0.005, p(3), p(4)/3]);

% ah{3,1}(3): right panel to (3,1)
  set(ah{3,1}(3),'Position',[p(1)+p(3)+0.005, p(2), p(3)/6, p(4)]);

setax(ah{1,1}(1))

% cg=[1 1 1]*.65;
  cg=[.8 .6 .4];

  lo={{'Color','k'},{'Color',cg}};

  if ~isempty(EE)
     n=size(EE,2);
     i=2:2:n; h=plot(i,EE(1:nE,i),lo{2}{:}); set(h(1),'Disp','odd'); hold on
     i=1:2:n; h=plot(i,EE(1:nE,i),lo{1}{:}); set(h(1),'Disp','even');
     ylim([0 2*Inrg.ops.ff(end)]);
  end

  rr=cat2(1,II(:,end).rr); rr(find(rr<=0))=nan;
  rr=-log(rr); nr=min(512,size(rr,2));

  legdisp('Location','NorthEast','Orientation','horizontal','eraset','pick',[2 1]);
  ylabel(['energy' char(10) 'flow diagram']);
  set(gca,'XTickLabel',[]);

setax(ah{1,1}(2))

  lo={{'Color','k'},{'Color',cg}};

  n=size(rr,1);
  i=1:2:n; h=plot(i,rr(i,1:nr),lo{1}{:}); set(h(1),'Disp','even'); hold on
  i=2:2:n; h=plot(i,rr(i,1:nr),lo{2}{:}); set(h(1),'Disp','odd');

  sms(2); xtight; ylim([-0.05 5]);
  legdisp('Location','NorthEast','Orientation','horizontal','eraset');

  ylabel(['entanglement' char(10) 'flow diagram']);
  set(gca,'XTickLabel',[]);

setax(ah{2,1}(end))

  h1=semilogy(e1); hold on
  set(h1(1),'Disp','e1: <eig(\rho_n)> using \chi_N'); % mv2front(h1(1));

  h2=semilogy(e2); % mv2back(h2);
  set(h2(1),'Disp','e2: <eig(\rho_n)> using \chi_E'); % mv2front(h2(1));

  blurl(h1,'cfac',0.85); % mv2front(h1(1));

  if ~isempty(yl), ylim(yl);
  else ytight(1.2); y=ylim; if y(2)<1E-4, y(2)=1E-4; ylim(y); end
  end
  set(ah{2,1}(end),'XLim',get(ah{1,1}(1),'XLim'),'YTick',logspace(-15,0,6));

  label('Wilson iteration n','\epsilon_n^D');
% mvlabel('y',0.02);

% leg21:bottom
  legdisp('Location','SE',... %'Orientation','horizontal',...
    'dx',[ 0.01  0.01],'dy',[-0.2],'ysc',0.6,'fs',10,'-detach');
  h=ymark(dw,'r--','-fg','istr',{'\epsilon^D',0.10,'VerticalAl','bottom'});

setax(ah{3,1}(1)); j=1; n=size(II,1); h=zeros(n,1);

  ex=zeros(n,1);
  for i=1:n, setax(ah{3,1}(1)); j=1;
% for j=1:size(II,2) % pretty much repeats the same data
     erd=II(i).erd; if isempty(erd) || i<k0+1, continue; end
     j=find(erd(:,2)>0); % rm(i)=min(erd(j,2));

   % extremal energies
   % ex(i)=max(erd(:,1));
     ex(i)=mean(erd(find(erd(:,1)>0.8*max(erd(:,1))),1));

   % plot eigenvalue without having degeneracy multiplied
     h(i)=semilogy(erd(j,1),erd(j,2),'.','Color',getcolor(i));
     hold on
  end

setax(ah{3,1}(1));

  k=find(h); h=h(k); ex=ex(k);

  if isempty(nx)
     ix=[ find(ex==min(ex)), find(ex==max(ex)) ];
   % sx={'n_{min}', 'n_{max}' };
     sx={'n_{max}', 'n_{min}' };
   % refers to max/min disc. weight // Wb,May11,14
  else
     ix=zeros(size(nx));
     for i=1:numel(nx), ix(i)=find(k==nx(i)); end
     sx={ 'n_1', 'n_2' };
  end

  om={'Color',[1 .3 .3],'LineW',2}; o={om{1:2},'Marker','o'}; % 'ro'
  i=ix(1);
    set(h(i),o{:},'Disp',sprintf('%s=%g',sx{1},k(i)));
    sms(h(i),3);

  om={'Color','k','LineW',2}; o={om{1:2},'Marker','d'}; % 'kd'
  i=ix(2);
    set(h(i),o{:},'Disp',sprintf('%s=%g (N=%g)',sx{2},k(i),Inrg.N));
    sms(h(i),3);

  mv2front(h(ix));
  hx=h; hx(ix)=[]; blurl(hx,'cfac',0.7);

  label(['rescaled energy E^{[n]}' bs], ['\rho^{[n;n_0]}' bs]);
% mvlabel('y',0.04);
  set(gca,'YTick',logspace(-20,0,5));

setax(ah{2,1}(1));

  n=size(II,1); lo={'Color',[1 1 1]*.65};
  i=2:2:n; h1=plot(i,(se(i,:)),lo{:});   set(h1(1),'Disp','even'); hold on
  i=1:2:n; h2=plot(i,(se(i,:)),'--',lo{:}); set(h2(1),'Disp','odd');

  h=[h1(end),h2(end)];
  ds='n_0\rightarrow \infty';
  set(h,'Color','b','LineW',2); set(h(1),'Disp',ds); mv2front(h);

% leg21:center
  legdisp('Location','NorthEast','Orientation','horizontal',...
    'dy',[0 0 -0.2]); %,'pick',[1 3 2] 
  y=ytight(1.1,'y1',1); ylabel('entropy S'); % mvlabel('y', 0.04)
  if y(2)<3,
       set(gca,'YTick',0:.5:y(2));
  else set(gca,'YTick',0:y(2)); end

setax(ah{2,1}(2));

  plot(Inrg.NK(:,1),'o-'); sms(4);
  ylabel('N_{keep}'); mvlabel('y',0.02);

  ytight(1.1,'y1',0);

  y=get(gca,'YTick'); if numel(y)<=2,
  set(gca,'YTick',linspace(0,y(end),3)); end

  h=ah{1,1}(1);
  set([ah{1,1}; ah{2,1}],'XLim',get(h,'XLim'),'TickLength',[0.005 0.005],...
    'XTick',get(h,'XTick'), 'XTickLabel',get(h,'XTickLabel'),'FontS',10);

  setax(ah{2,1}(end));
  set(gca,'XTickLabel',get(gca,'XTick'));

setax(ah{3,1}(1));

  ERD=cat(1,II.erd);
% xb=linspace(0,ceil(max(ERD(:,1))),32);
  n=ceil(22/sqrt(Lambda)); % n=16 for Lambda=2, n=8 for Lambda=8
  xb=linspace(0,Etr,16);

  wbhist([],[],xb); n=size(II,1); Np=0;
  for i=1:n
     erd=II(i).erd; if isempty(erd), continue; end
   % NB! do not include degeneracy in polynomial fit below!
     wbhist(erd(:,1),erd(:,2)); Np=Np+1;
  end

  [x,y]=wbhist; y=(y/Np)./diff2(x,'len');

  s='\rho(E)'; dwc=nan;

  if 1
  % i=find(y>1E-12); % q=log(y(i)); i=find(y>exp(mean(q)+std(q)));
  % q=max(x); i=find(x>0 & x<0.80*q);
    i=2:length(x)-1;

    i=i(find(y(i))); % plot(x,y,'k');
    p=polyfit(x(i),log(abs(y(i))),1); xx=xlim; xx=linspace(xx(1),xx(2));

    h1=plot(x(i),abs(y(i)),'r','Disp','cumulative weight');
    h2=plot(xx,exp(polyval(p,xx)),'k');
    mv2back(h2); blurl([h2]);

  % ensure that rho(E) is normalized function for E in [0,infty]
    q=p; q(2)=log(abs(q(1)));
    h=plot(xx, exp(polyval(q,xx)),'-.','Color',[.3 .3 .3],...
      'Disp',sprintf('%s \\sim e^{%.4g{\\cdot}E}',s,p(1)));

    dwc=exp(p(1)*Etr);


  else
    i=find(ERD(:,2)>1E-12);
    i=find(ERD(:,2)>exp(mean(log(ERD(i,2)))));

    p=polyfit(ERD(i,1),log(ERD(i,2)),1);

    xx=xlim; xx=linspace(xx(1),xx(2));
    h=plot(xx,exp(polyval(p,xx)),'k','Disp',...
      sprintf('%s \\approx %.3g e^{%+.4g{\\cdot}E}',s,exp(p(2)),p(1)));
  end
  header('fleft','-append',' \rho(E) fit: %.4g',p(1))

  Iout.p=p;

% q=3;
% h=plot(xx,exp(-q*sqrt(Lambda)*xx),'m--','Disp',...
%   sprintf('%s \\sim e^{-{%g}\\Lambda^{1/2}E}',s,q));

  if ~isempty(xl31)
     if numel(xl31)==2, xlim(xl31); else axis(xl31); end
  else xtight(1.05); %'x1',-0.25);
  end

  if ~isempty(yl31)
     if numel(yl31)==2, ylim(yl31); else axis(yl31); end
  else ytight(1.1,'y1>',1E-15); end

  o={'-detach','fs',10};
  if nochi
     legdisp('pick',[3 4],o{:},'dx',[0.02 0.005],'dy',[0 0.05]); y=0.65;
  else
     ymark(dw,'r--','-fix','Color',[1 .4 .2],'-fg');
     xmark(Etr,'r--');

   % leg31
     legdisp('pick',[3 4],'dy',[0 0.10],o{:},'ysc',0.7),'-detach';
  end

  if ~isfield(Inrg,'Itr') || isempty(Inrg.Itr) || all(Inrg.Itr.Etrunc<=0)
     s=sprintf('\\Lambda=%g, N_K=%g\nn_0=%g',Inrg.Lambda,Inrg.Nkeep,n0_);
  else
     q=Inrg.Itr.Etrunc(1); % Etrunc(2) is relEtrunc if any
     e=norm(diff(Inrg.ops.ff(end-3:end,:)))/abs(Inrg.ops.ff(end));
     if e>1E-2, error('Wb:ERR','\n   ERR couplings still vary'); end
   % Efac=1/abs(Inrg.ops.ff(end)); q=q*Efac;
     s=sprintf('\\Lambda=%g,E_K=%.2g\nn_0=%g',Inrg.Lambda,q,n0_);
  end
  postext({'E',[0.02 0.12]},s,{'FontS',10});

% leg31
  legdisp('Location','SouthWest','pick',[2 1],'fs',10,'eraset',...
    'dx',[-0.02 -0.005],'dy',-0.05,'xsc',{0.8 [0 1.5]},'-detach','ysc',0.8);

setax(ah{3,1}(2)); % top panel to (3,1)

  xb=linspace(0,ceil(max(ERD(:,1))),32);
  wbhist([],[],xb); n=size(II,1);

  for i=1:n
% for j=1:size(II,2) % pretty much repeats the same data
     erd=II(i).erd; if isempty(erd), continue; end
     wbhist(erd(:,1),[]);
% end
  end

  [ee,nn,xb]=wbhist; nn=nn/Np;

  plot(ee,nn,'o-'); sms(3); hold on
  ylabel('\nu(E)');

  ytight(1.2,'y1',0);
  set(gca,'XLim',get(ah{3,1}(1),'XLim'),...
    'XTick',get(ah{3,1}(2),'XTick'),'XTickLabel',[]);

% header('SE',sprintf('%s//N_{tot}=%g',hostname,sum(nn)))

setax(ah{3,1}(3)); % right panel to (3,1)

  wbhist([],[],linspace(-60,0,48)); n=size(II,1);

  for i=1:n
% for j=1:size(II,2) % pretty much repeats the same data
     erd=II(i).erd; if isempty(erd), continue; end
     j=find(erd(:,2)>0);
     wbhist(log(erd(j,2)),[]);
% end
  end
  [logr,nn,xb]=wbhist; nn=nn/Np;

  semilogy(nn,exp(logr),'o-'); hold on
  sms(3); % scalecolor('-flip','alpha',2);

  if ~nochi
     h=ymark(dw,'r--','-fg'); x=xlim;

     if smooth>0
          s=sprintf('\n    [using smooth=%g]',smooth);
     else s=''; end

   % NB! adjust vertical position of arrow to actually match dw (roughly)
     text(x(2)+0.3*diff(x),2*dw, ['\rightarrow  ' sprintf(... 
     '\\epsilon^D_{\\chi=%g} = %s%s',chi,num2tex(dw,'%.2E'),s)], ...
     'Color','r','VerticalAl','top');

     h=ymark(dwc,'k--');
     text(x(2)+0.3*diff(x),2*dwc, ['\rightarrow  ' sprintf(... 
     '\\epsilon^D_{\\chi=%g} = %s%s',chi,num2tex(dwc,'%.2E'),s)], ...
     'VerticalAl','top');
  end

  xtight(1.2,'x1',0); 
  xlabel('\nu(\rho)');

  set(gca,'YLim',get(ah{3,1}(1),'YLim'),'YTick',...
  get(ah{3,1}(1),'YTick'),'YTickLabel',[]); %,'YAxisLocation','right');

  setax(ah{2,1}(end,1)); % mv2front

end

% -------------------------------------------------------------------- %
% -------------------------------------------------------------------- %

function [dd,neg]=seteps2nan(dd,eps)

   if nargin<2, eps=0; end

   if nargout>1
      i=find(dd<0); neg=numel(i); if neg
        neg=neg/numel(find(~isnan(dd))); % return rel. percentage
        dd(i)=abs(dd(i));
     end
   end

   i=find(dd<=eps); e=norm(dd(i));
   if e>1E-12
      s=sprintf('WRN skipping small weight e=%g !??',e);
      if e<1E-6, wblog('WRN',s'); else
      error('Wb:ERR','\n   ERR %s',s); end
   end

   dd(i)=nan;

end

% -------------------------------------------------------------------- %
% -------------------------------------------------------------------- %

