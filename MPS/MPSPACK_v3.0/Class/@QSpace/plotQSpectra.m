function Iq=plotQSpectra(H,varargin)
% function Iq=plotQSpectra([H|R] [,opts][,SOP])
%
%    graphical analysis of eigenspectrum of Hamiltonian
%    (or entanglement spectrum of density matrix; requires '-ES' flag)
%    separating the data into its symmetry subspaces.
%
% Options
%
%   '-ES'     entanglement spectrum (expects Rho as input)
%   '-x'      show state space dimensions (rather than multiplet dimensions)
%
%   'yl',...  ylim range
%   'y2',...  y-range within which degeneracies are indicated
%   'ws',...  which symmetries to pick (e.g. [1 2] for the first two)
%
%   SOP       more verbose description of symmetry labels, which also
%             appears in the title of the respective panels.
%             SOP my simply be specified through IS.SOP, i.e. as returned
%             by getLocalSpace() in the NRG setup; all that is used, however,
%             are the strings SOP(:).info, hence, alternatively, SOP may
%             only contain {SOP.info}. keywords looked for are 'spin' or
%             'charge' which effects the xlabels.
%
% Adapted from getrhoESpectra.m
% Wb,Apr11,12

   getopt('INIT',varargin);
      Rflag= getopt('-ES');
      scale= getopt('-sc');
      e0   = getopt('e0',[]);

      if     getopt('-x'), vflag=1;
      elseif getopt('-v'), vflag=2;
      elseif getopt('-V'), vflag=3;
      else vflag=0; end

      yl = getopt('yl',[]);
      y2 = getopt('y2',[]);
      is = getopt('ws',[]);
      mat= getopt('mat','');
      istr=getopt('istr',-1);
    % IST= getopt('istr',I.IST);
      eps= getopt('eps',0.1);
      ah = getopt('ah',[]);
      odeg=getopt('odeg',{});
 % getopt('check_error');
   SOP=getopt('get_last',[]);

   if isempty(H)
      error('Wb:ERR','\n   ERR invalid usage (empty H)');
   end

 % choose rank-1 symmetry sectors by default.
 % if there are no more than one abelian symmetry sector,
 % look out for rank-2 symmetries
   d=getsym(H,'-d');
   if isempty(is)
      i=find(d==1); if numel(i)>1, is=i(1:2);
      elseif numel(d)>1
        is=[i, find(d==2)]; if numel(is)>2, is=is(1:2); end
      else is=1;
      end
   end

   ns=numel(is);
   for i=1:ns
      ws(i)=getsym(H,'-I',is(i)); % disp(ws(i))
   end

   q=isdiag(H);
   if ~q
      [ee,Ie]=eigQS(H); R_=H; q=2; ee=ee(:,1);
      if ~isreal(ee)
         error('Wb:ERR','\n   ERR got complex eigenvalues !??'); end
      H=QSpace(Ie.EK);
   else
      if q==1
         for i=1:numel(H.data), H.data{i}=diag(H.data{i}).'; end
         q=2;
      elseif q~=2, error('Wb:ERR',['\n   ' ...
         'ERR invalid input QSpace (not an operator!? [%g])'],q);
      end
      ee=sort(cat(2,H.data{:}));
   end

   if Rflag
      if any([abs(sum(ee))-1, -sum(ee(find(ee<0)))]>1E-12)
         error('Wb:ERR','\n   ERR got invalid density matrix'); end
      ee(find(ee<=0))=[]; ee=sort(-log10(ee)); e0=min(ee);
      for i=1:numel(H.data), dd=H.data{i}; 
         dd(find(dd<=0))=[];
         H.data{i} = sort(-log10(dd))-e0;
      end
      ylb='entanglement spetrum (using log_{10})';
   else
      if isempty(e0), e0=min(ee); end
      if abs(e0)>1E-12
         for i=1:numel(H.data)
            H.data{i}=H.data{i}-e0; 
         end
      end
      ylb='energy';
   end

   if nargout
      Iq=add2struct('-',H,ee,e0,Rflag);
   end

  if isempty(ah)
     s=mfilename; if Rflag, s=[s ':R']; else s=[s ':E']; end
     ah=smaxis(1,ns,'tag',s,'zoom',[0.95 0.8],'dx',0.02);
     if numel(ah)>1
        p=get(gcf,'Pos'); q=350*ns;
           p([1 3])=[sum(p([1 3]))-q, q]; q=400;
           p([2 4])=[sum(p([2 4]))-q, q];
        set(gcf,'Pos',p);
     end
     hflag=1;
  else hflag=0;
     if numel(ah)<ns, error('Wb:ERR',['\n   ERR ' ...
       'insufficient number of axis handle (%d/%d)'],numel(ah),ns); 
     end
  end

  if hflag || vflag
     s={};
     if e0~=0
        s{end+1}=sprintf('shift=%.4g',e0);
     end
     if isfield(H.info,'qtype')
        s{end+1}=sprintf('''%s''',getsym(H));
     end
     if isfield(H.info,'itags') && ~isempty(H.info.itags)
      % s{end+1}=[ itags2str(H.info.itags) ];
        s{end+1}=[ '''' H.info.itags{1} '''' ];
     end
     s=strhcat(s,'-s',', ');

     if Rflag
          s=['%M::entanglement spectra: ' s];
     else s=['%M::eigenspectra: ' s]; end


     if isequal(istr,-1)
        if ~isempty(mat), s=sprintf('%s (%s)',s,mat); end
     elseif ~isempty(istr)
        s=sprintf('%s [%s]',s,istr);
     end

     header(s); addt2fig Wb

     if vflag==1 && isfield(H.info,'cgr')
        setax(header), postext([0.01 0.93],['NB! showing ' ... 
       'full degeneracies of state spaces (rather than multiplets)'],...
       {'FontSize',10,'Color',[.6 .2 .1]});
     end
  end

  r2=0.40; lo={'Color',[.7 .7 .7]};

  if isempty(y2)
     if Rflag, y2=0.4*max(ee); else y2=0.6*max(ee); end
  end

  if isempty(yl)
     yl=max(ee)-min(ee);
   % yl=e0+yl*[-0.02, 0.8];
     yl=   yl*[-0.02, 0.8]; % Wb,Mar30,16
     if Rflag, yl(1)=-0.02*diff(yl); end
  elseif numel(yl)==1
     q=max(ee)-e0; % min(ee);
     yl=[ min(ee)-0.02*q, yl ];
  end
  if ~isempty(yl) && y2>yl(2), y2=yl(2); end

  dE=[];
  if isempty(odeg), odeg={ 'dx',[0.18 0.35],'n',256,'y2',y2 }; end

  for k=1:ns, setax(ah(k)); fs=[];
     [q,EE,qq,dz,DZ]=getQSpectra(H,ws(k)); nq=size(q,1);
     if nargout
        Iq.q{k}=q; Iq.EE{k}=EE; Iq.qq{k}=qq; Iq.dz{k}=dz;
     end

     if size(q,2)>1
      % xq=max(q,[],2); xq=0.5*xq.*(xq+1);
      % xq=q(:,1).^2 + 2*q(:,1).*q(:,2) +q(:,2).^2;
      % i=find(q(:,1)>q(:,2));
      % xq(i)=-xq(i);
      % ------------------------------------------------------------ %
      % label from smallest to largest multiplet
      % where scalar representation is at xq=0 // Wb,Sep07,15
      % in case of degeneracy, sort with respect to q, while taking
      % into account dual representations => sort(q,2) // Wb,Sep26,15
      % ------------------------------------------------------------ %
      % [~,i1,i2]=unique(dz); xq=0:length(i1)-1; xq=xq(i2);
        [~,is]=sortrows([dz,sort(q,2)]); xq(is)=0:length(is)-1;  ql=q;
      % x=q-fliplr(q);
      % distinguish between dual representation
      % for i=1:size(x,1)
      %    j=find(x(i,:),1); if ~isempty(j) && x(i,j)<0, xq(i)=-xq(i); end
      % end
     else xq=q; end

     for i=1:nq, x=xq(i); ee=EE{i}; if isempty(ee), continue; end
        xx=repmat([x-r2;x+r2;nan],1,numel(ee)); xx(3,:)=nan;
        yy=repmat(ee,3,1);
        plot(xx(:),yy(:),lo{:}); hold on

        if vflag==1, osd={'dz',DZ{i}};
        elseif vflag==2, osd={'dz',DZ{i},'-v'};
        elseif vflag==3, osd={'dz',DZ{i},'qq',qq{i}};
        else osd={}; end

        show_deg(ee,gca,osd{:},'xt',x,odeg{:});
      % h=ymark(ee); set(h,'Color',[.7 .7 .7],'LineW',1,'LineSt',':');
     end

     if ~isempty(yl), ylim(yl); end
     xtight(1.1,'view')
     h=xmark(0);h(2)=ymark(0); set(h,'Color',[1 .9 .8]);
     set(h(1),'LineW',12); set(h(2),'LineW',2);

     if isset('ql')
        n=size(ql,1); s=repmat({''},1,n);
        for i=1:n
           l=['(' sprintf('%g',ql(i,:)) ')'];
           if isempty(s{i}), s{i}=l; else s{i}=[s{i} 10 l]; end
        end
        [xq,is]=sort(xq);
        set(gca,'XTick',xq,'XTickLabel',s(is));
        fs=8;
     end

     if ~isempty(SOP)
        if iscell(SOP)
             s=SOP{ws(k).is};
        else s=SOP(ws(k).is).info; end

        title(s); ic=[];
        switch ws(k).sym
           case 'A'
              if ~isempty(findstr(s,'charge')), xlabel('Q'); ic=k;
              elseif ~isempty(findstr(s,'spin')), xlabel('2S_z'); end
           case 'SU2'
              if ~isempty(findstr(s,'charge')), xlabel('2Q');
              elseif ~isempty(findstr(s,'spin')), xlabel('2S'); end
           otherwise
        end
        if ~isempty(ic)
         % find charge excitation energy for same ground state Q
           e0=min(cat(2,EE{:})); k0=k;
           for i0=1:nq
              if ~isempty(find(EE{i0}<=e0)), break; end
           end
           de=diff(EE{i0}); dE=de(find(de>0.05,1));
        end
     else
        s=ws(k).sym;
        if isequal(s,'A'), title('abelian');
        elseif isequal(s,'SU2'), title('SU(2)'); xlabel('2S');
        elseif ~isempty(regexp(s,'SU\d'))
           n=str2num(s(3:end));
           title(sprintf('SU(%g)',n));
           xlabel('multiplet label (sorted by dimension)');
        end
     end

     if k==1, ylabel(ylb); end
     if isset('fs'), set(gca,'FontSize',fs); end
  end

  if ~isempty(dE) && scale
     for k=1:ns, setax(ah(k));
        scaley(dE,'ulabel','{\delta}E');
        if k>1, ylabel(''); end
     end
     setax(ah(k0)); grid on
  end

end

% -------------------------------------------------------------------- %
% -------------------------------------------------------------------- %

