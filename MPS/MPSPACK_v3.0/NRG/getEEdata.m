function [EE,hh,iQ]=getEEdata(HK,varargin)
% function [EE,hh,iQ]=getEEdata(HK [,opts])
%
%    Auxilliary routine for plotting NRG like energy flow diagrams:
%    combine energies from give set of QSpaces while also grouping
%    symmetry spaces(!).
%
% Options
%
%   'EK',..  keep only data up to energy EK
%   'E0',..  E0 data from NRG data to determine even/odd sector for ground state
%   '-red'   reduce to unique data
%   '-rel'   relative, ie subtract lowest energy
%   '-ES'    plot entanglement spectra (assuming got density matrices as input)
%
%   'ah',..  specify axis handles to plot into (automatically sets '-p')
%   'dk',..  specify interval length and offset: 'dk',[dk [,k1]]
%   '-p'     plot result
%   '-1'     do not distinguish between even or odd (plots into single axis)
%            equivalent to: 'dk',[1]
%   '-eo1'   plot into single axis, yet do distinguish between even and odd
%            equivalent to: 'dk',[2],'ah',<single axis handle>
%   'yl',..  specifies ylim (data fully outside this window is skipped)
%   '-ldisp' set 'Disp' for legend for each symmetry sector
%
%   'Qmap',{qm1,qm2,...}
%        allow prior map on Q-labels in terms of linear superposition
%        or additive terms (possibly as string function of iteration index 'i')
%        where the order may change (hence allow for general sequence).
%    Ex: switching from particle (n1,n2) to
%        (charge,spin) = (n1+n2-1, n2-n1) <=> Qmap = { [ 1, 1; -1 1 ], '[-i;0]' }
%        1st: rotate to total charge (N) and total spin (Sz)
%        followed by 2nd: shift to relative to half-filling
%
% Wb,May29,11

  getopt('init',varargin);
     Ekeep =getopt('EK',Inf);
     E0    =getopt('E0',[]);
     qsel  =getopt('qsel',[]);
     pflag =getopt('-p');
     dk    =getopt('dk',[]); % Wb,Sep07,15
     kQ    =getopt('kQ',[]); % Wb,Sep26,15
     cdata =getopt('c',{});
     redfl =getopt('-red');
     rflag =getopt('-rel');
     ESflag=getopt('-ES');   % Wb,Dec13,16
     eps   =getopt('eps',1E-2);
     dx    =getopt('dx',0);
     yl    =getopt('yl',[]);
     k2x   =getopt('k2x',{}); % { x(k), 'xlabel' }
     Qmap  =getopt('Qmap',{});

     ah=getopt('ah',[]);
     if isempty(ah) && getopt('-gca'), ah=gca; end
     ldisp =getopt('-ldisp');

     if getopt('-1')
      % single plot (e.g. for entanglement spectra)
        nax=1; dk=1; if isempty(ah), ah=gca; end
     elseif getopt('-eo1')
        nax=1; dk=2; if isempty(ah), ah=gca; end
     elseif ~isempty(dk)
        if numel(dk)==1, nax=dk;
        else nax=1; end
     elseif ~isempty(ah)
        dk=numel(ah); nax=dk;
     else nax=0; end

     if nax && ~pflag, pflag=2;
     elseif ~nax && pflag, nax=1; dk=1; end

  getopt('check_error');

  if ~isempty(HK(1).data)
     isd=isdiag(QSpace(HK(1)),'-d');
     if isd>2
        wblog('NB!','using H''');
        HK=HK';
     end
  else isd=-1; end

  HK=struct(HK);
  L=numel(HK); if L<3, error('Wb:ERR',...
    '\n   ERR invalid HK (got only QSpace vector of length %g)',L); 
  end

  if ~isempty(k2x)
     if ~iscell(k2x) || numel(k2x)<2 || ~ischar(k2x{2})
        error('Wb:ERR','\n   ERR invalid usage (k2x)');
     end

     s=k2x{2}; if numel(k2x)>2, k1L=k2x{3}; else k1L=[1 L]; end
     k2x={k2x{1}, k1L, s};
  end

  if ~isempty(Qmap)
     if ~iscell(Qmap), Qmap, error('Wb:ERR','\n   ERR invalid Qmap'); end
     qdim=size(HK(1).Q{1},2);
     for k=1:numel(Qmap), 
      % allow to include mathematical expression in terms of i below
        [qm,qms]=test_qm_string(Qmap{k}); s=size(qm);

        if isequal(s,[qdim,qdim])
           for i=1:L, Q=HK(i).Q;
              if ~isempty(Q), m=numel(Q); else continue; end
              if ~isempty(qms), eval(['qm=' qms ]); end
              for j=1:m, Q{j}=Q{j}*qm; end
              HK(i).Q=Q;
           end
        elseif isequal(s,[1,qdim])
           for i=1:L, Q=HK(i).Q;
              if ~isempty(Q), m=numel(Q); else continue; end
              if ~isempty(qms), eval(['qm=' qms ]); end
              for j=1:m, Q{j}=Q{j} + repmat(qm,size(Q{j},1),1); end
              HK(i).Q=Q;
           end
        else Qmap, error('Wb:ERR','\n   ERR invalid Qmap'); end
     end
  end

  if isd==0 % && ~nax
   % HK(1) may simply be the non-diagonalized H0
   % Wb,Apr13,13 // tags: KEEP_ITER0_AH
   % Wb,Nov06,15 :: commented && ~nax above: *always* diagonalize!
   % e.g. required by: if nargout => EE=cat(2,EE{:}); below
     wblog('NB!','diagonalizing H0');
     [ex,I]=eigQS(QSpace(HK(1)));
     HK(1)=QSpace(I.EK) - min(ex(:,1));
  end

  if ESflag
     for k=1:numel(HK), Hk=HK(k).data;
        for j=1:numel(Hk)
           Hk{j}=fliplr(-log10(Hk{j}(find(Hk{j}>1E-14))));
        end
        HK(k).data=Hk;
     end
  end

  QQ=cell(1,L); DD=cell(1,L); DC=cell(1,L);
  for i=1:L, H=HK(i); if isempty(H.Q), continue; end
   % truncate energy data if Ekeep<Inf is specified
     if Ekeep<Inf, n=numel(H.data);
       for j=1:n, h=H.data{j};
          H.data{j}=h(find(h<=Ekeep));
       end
       HK(i)=H;
     end
     [QQ{i},DD{i},DC{i}]=getQDimQS(H,2);
  end

  Q=uniquerows(cat(1,QQ{:})); nQ=size(Q,1); EE=cell(nQ,L);
  D=zeros(1,nQ); d=[]; % zeros(nQ,size(DC{1},2));

  for i=1:L, if isempty(QQ{i}), continue; end
     [i1,i2,I]=matchIndex(Q,QQ{i});
     if ~isempty(I.ix2), error('Wb:ERR','\n   ERR missing Q-data !??'); end
     D(i1)=max(D(i1),DD{i}(i2));
     d(i1,:)=DC{i}(i2,:);
     EE(i1,i)=HK(i).data(i2);
  end

  EE_=EE;

% combine data for each symmetry sector across entire system
% into block matrix (thus keep site resolution)
  for i=1:nQ, EE{i}=cat2(1,EE{i,:},{nan}); end
  EE=EE(:,1);

% sort symmetry sectors with respect to energy
  if isempty(kQ), kQ=max(ceil(L/2),L-3); end

  if ~isempty(dk)
       k=ceil(dk/2); k=kQ+(-k:+k);
  else k=kQ:kQ+1; end
  q=max(k)-L; if q>0, k=k-q; end

  E0=nan(nQ,1);
  for i=1:nQ, if ~isempty(EE{i})
     q=EE{i}(k,:); E0(i)=min(q(:));
  end, end

% [D,is]=sort(D,'descend');
  [E0,is]=sort(E0); D=D(is);
     ix=find(D==0); if ~isempty(ix), D(ix)=[]; is(ix)=[]; nQ=numel(D); end
  EE=EE(is,1); Q=Q(is,:); d=d(is,:);

  if ~isempty(qsel)
     if size(qsel,2)~=size(Q,2)
        error('Wb:ERR','\n   ERR invalid qsel (size mismatch)'); 
     end
     qsel=uniquerows(qsel);
     [ia,ib,Im]=matchIndex(Q,qsel);
     if isempty(ia), Q, qsel
        error('Wb:ERR','\n   ERR invalid qsel (non-matching Q)');
     end
     wblog(' * ','using qsel (%g)',size(qsel,1));
     EE=EE(ia); Q=Q(ia,:); nQ=size(Q,1);
     D=D(ia); d=d(ia,:);
  end

  if redfl, nx=0; n1=sum(D);
     for i=1:nQ, dd=diff(EE{i},[],2); dd(find(isnan(dd)))=0;
         j=find(sum(dd.^2,1)<eps);
         if ~isempty(j), nx=nx+numel(j);
             EE{i}(:,j+1)=[]; D(i)=size(EE{i},2);
         end
     end
     if nx, n2=sum(D);
        wblog(1,'<i> removing %g degenerate records (%g, E<%g: %.3g%%)',...
        nx,eps,Ekeep,100*nx/n1)
     end
  end

  if rflag
     ee=min(cat(2,EE{:}),[],2);
     for i=1:nQ, EE{i}=EE{i}-repmat(ee,1,size(EE{i},2)); end
  end

  if ~pflag
   % return [EE,Q [,iQ]] // Wb,Feb27,17
     if nargout>1, hh=Q; end
     if nargout>2
        n=numel(EE); iQ=cell(1,n);
        for i=1:n, iQ{i}=repmat(i,1,size(EE{i},2)); end
        iQ=[iQ{:}];
     end
     if nargout~=2, EE=cat(2,EE{:}); end
     return
  end

  nah=numel(ah);

  if nax>1
     if nah==1
      % ah=repmat(ah,1,nax); // keep unique! // Wb,May17,17
     elseif isempty(ah)
        ah=smaxis(nax,1,'tag',mfilename);
        header('%M'); addt2fig Wb
     elseif ~all(isaxis(ah)) || numel(ah)<nax
        error('Wb:ERR','\n   ERR getEEdata: got invalid axes handles');
     end
  end

  nc=numel(cdata); % number of colors to use

  ns=dk(1);
  if numel(dk)==2, jset=dk(2); else jset=1:ns; end
  nj=numel(jset);

% NB! symmetries may appear in even XOR odd symmetry sector
  mm=ones(nQ,ns);
  for i=1:nQ, dd=EE{i}; n=size(dd,1); 
     for j=1:ns
       if all(isnan(reshape(dd(j:ns:end,:),[],1))), mm(i,j)=0; end
     end
  end

  for j=1:ns, i=find(mm(:,j));
     mm(i,j)=1:numel(i); % color index // sort
  end

  co=get(ah(1),'ColorOrder');
  if isequal(co,get(0,'defaultAxesColorOrder'))
     co(end,:)=[1 1 1]*0.6; % default: [ 0.25 0.25 0.25 ]
   % co=[0 0 0; co]; % include black as 1st color
   % co=[co; 0 0 0]; % include black as last color
   % NB! stick with default color order to start with! // Wb,Sep26,15
     co=[co
        0 0 0 % black     % see also MLIB/getcolor.m
        0.88  0.68  0.3   % brown
        0.5   0.75  1     % light blue
        0.88  0.68  0.34  % orange.in
     ];
  end

  setuser(gcf,'cmap',co);
  setuser(gcf,'Q',Q);
  setuser(gcf,'D',D);
  setuser(gcf,'d',d);

  Iq=struct('sym',[],'Q',[],'i',[],'j',[],'o',[]);
  if ~isempty(HK(1).info)
       Iq.sym=HK(1).info.qtype;
  else Iq.sym='A*'; end

  hh=cell(nQ,ns); bflag0=(nax==1 && numel(dk)<2);

  iQ=Q; ih=size(Q,2)+1;

  for iq=nQ:-1:1, dd=EE{iq}; n=size(dd,1); 
   % plot largest (ground state space) last
     if nc, c=cdata{mod(iq-1,nc)+1};
   % elseif nax>1, if mod(iq,2), c='k'; else c='r'; end
     else c=[]; end

     Iq.Q=Q(iq,:); Iq.iq=[iq nQ];

     for j1=1:nj, j=jset(j1); % 1:ns
        if mm(iq,j)
           if nax>nah % Wb,May17,17
              l=mod(j1-1,nah)+1; bl=(nj-j1)/nj;
           else l=j1; bl=0; end
           setax(ah(l));

           if isempty(c)
              % o=getlopts(mm(iq,j),co);
                o=getlopts(iq,co); % be consistent across panels // Wb,Sep26,15
           else o={'Color',c}; end
           if bflag0 && j==2, bflag=1; % blurl
              o{2}=1-0.3*[1-o{2}];
           else bflag=0; end

           ij=j:ns:n; x=ij+dx; y=dd(ij,:); [nx,ny]=size(y); hl=[]; ym=inf(ny,1);
           for l=ny:-1:1
              i=find(~isnan(y(:,l)));
              if numel(i)>min(0.9*nx,nx-3) % Wb,May29,17
               % if there are only certain individual iterations with nan
                 i=1:nx;
              end

              if ~isempty(i) && ~isempty(yl)
                 ir=find(y(i,l)>=yl(1) & y(i,l)<=yl(2)); % data within y-range
              else ir=1; end
              if ~isempty(i) && ~isempty(ir)
                 hl(end+1)=plot(x(i),y(i,l),o{:}); hold on
                 ym(numel(hl))=mean(y(i,l));
              end
           end
           hh{iq,j}=fliplr(hl);

           if ~isempty(hl), l=find(ym==min(ym),1); % ym(l)
              if ~isempty(l) % Wb,Jul19,17
                 iQ(iq,ih)=hl(l);
                 if ldisp && ~bl
                    set(hl(l),'Disp',sprintf('(%s)',sprintf('%g',Iq.Q)));
                 end
              end
           end

           if bl, blurl(hl,'cfac',bl,'LineWidth',1); end

           Iq.j=[j ns]; Iq.o=o;
           set(hh{iq,j},'UserData',Iq,'tag',sprintf('[%g]',Iq.Q));

           if bflag, mv2back(hh{iq,j}); end
        end
     end
  end

  for j=1:ns, hh{1,j}=cat(2,hh{:,j}); end
  hh=hh(1,:);

  set(ah,'XLim',[1 L]);
  if ~isempty(k2x)
     for i=1:numel(ah), setax(ah(i))
        map_xdata_1(k2x{:});
     end
  end

  if isempty(yl), ytight(1.1);
   % e2=[Inf,-Inf];
   % for i=1:numel(EE)
   %    e2(1)=min(e2(1),min(EE{i}(:)));
   %    e2(2)=max(e2(2),max(EE{i}(:)));
   % end
  else set(ah,'YLim',yl);
  end

  if nargout && nargout<3, EE=cat(2,EE{:}); end

end

% -------------------------------------------------------------------- %
% -------------------------------------------------------------------- %

function [qm,qms]=test_qm_string(qm)

  if ~ischar(qm), qm=qm'; qms=''; return; end

% allow to include mathematical expression in terms of index i only
  i=0; qms=[qm ''';'];
  eval(['qm=' qms]);

end

% -------------------------------------------------------------------- %
% -------------------------------------------------------------------- %

function map_xdata_1(k2x,k1L,s)

  if numel(k2x)==2
   % interpreted has mapping k=[1 L] to k2x([1 2])
     map_xdata(gca,'lin',k1L-1,k2x); % Wb,Sep13,15
  else
   % assuming that the current x-data is plain integer values (k)
     hh=findall(gca,'type','line')';
     for h=hh
        x=get(h,'XData'); j=find(~isnan(x));
        if ~isempty(j), x(j)=k2x(x(j)); set(h,'XData',x); end
     end
  end

  xlabel(['bond index k  \rightarrow  ' s]);
  xtight

end

% -------------------------------------------------------------------- %
% -------------------------------------------------------------------- %

