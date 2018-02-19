
% -------------------------------------------------------------------- %
% -------------------------------------------------------------------- %

  if ~exist('param','var'), global param; end
% e.g. when loaded from NRG_info file
  if isempty(param)
     param=add2struct('-',Lambda,'epsd?','Gamma?','Nkeep?');
  end

% if exist('eps')~=1 && exist('Inrg') && isfield(Inrg,'TN')
% eps=50*Inrg.TN; end

  if isset('noraw'), osmo={'raw'}; else osmo={'RAW'}; end
  setopts(osmo,...
    'afac?','sigma?','emin?','emax?','nlog?','reps?','eps?','skip?','keepA0?',...
    'xli?','yli?','adisp?');

  if ~exist('Idma','var')
     Idma.rho=nan(1,2);
     Idma.ver='';
     Idma.T=nan;
  end

  [ox,ax,ah,Idma.smo] = rsmoothSpec(om, a0, Idma, osmo{:});

  if exist('ac_cmd','var') && ischar(ac_cmd) && ~isempty(ac_cmd)
   % eg. ac_cmd='[ac,Ic]=getGC_JJH(ox,ax(:,1),-JH*ax(:,2),''-q'');'
     try
        e=[]; eval(ac_cmd);
     catch e % e=lasterror; 
        wblog(['ERR\N\N  ' ...
          'ERR failed to evaluate ac_cmd:\N  %s\N\N  msg: \\'],ac_cmd);
        disp(e.message); fprintf(1,'\n');
     end
     if isempty(e)
        if exist('ac','var') && size(ac,1)==length(ox) && isreal(ac)
           lo={'Color',[0 0 .4],'tag','Aimp'};
           for t={'Spec01','Spec01i','Spec02'}
              h=findall(gcf,'type','axes','tag',t{:}); if isempty(h), continue; end
              set(gcf,'CurrentAxes',h(1)); hold on
              h=addplot(ox,ac*Idma.smo.afac{1},lo{:});
              set(h(1),'Disp','A^{imp}(\omega)');
           end
        else 
           wblog('WRN','invalid ac_cmd (define variable ''ac'')');
        end

      % lo={'Color',getcolor('y2')};
        lo={'Color',[1 .6 .2],'tag','selfenergy'};
        % [1 .5 0]; 'm'=[1 0 1];

        if exist('Ic','var') && (exist('U','var') || exist('JH','var'))
        if isfield(Ic,'US')
           for t={'Spec01','Spec01i','Spec02'}
              h=findall(gcf,'type','axes','tag',t{:}); if isempty(h), continue; end
              set(gcf,'CurrentAxes',h(1)); hold on
              h=addplot(ox,-(1/(U*pi))*imag(Ic.US),lo{:});
              set(h(1),'Disp','\Sigma(\omega)');
           end
        elseif isfield(Ic,'SE')
           u=[]; if exist('U', 'var'), u(end+1)=U; end
                 if exist('JH','var'), u(end+1)=JH; end
           if ~isempty(u), u=max(u); else u=1; end

           for t={'Spec01','Spec01i','Spec02'}
              h=findall(gcf,'type','axes','tag',t{:}); if isempty(h), continue; end
              set(gcf,'CurrentAxes',h(1)); hold on
              h=addplot(ox,-(1/(u*pi))*imag(Ic.SE),lo{:});
              set(h(1),'Disp','\Sigma(\omega)');
           end
        else 
           wblog('WRN','invalid ac_cmd (define self-energy in ''Ic.US'')');
        end
        end
     end
  end

  nrg_header; TK=TKondo;

  s={sprintf('%s :: T=%.3g', Idma.ver, Idma.T)};
  if exist('B','var') && B~=0
     if TK~=0
          s{end+1}=sprintf(', B=%.3gT_K', param.B/TK);
     else s{end+1}=sprintf(', B=%.3g, T_K=%g', param.B, TK); end
  end

  header(cat(2,s{:}));
  addfinfo

% --------------------------------------------------------------------- %
% --------------------------------------------------------------------- %
% inset upper right
  setax(ah(3));

% inset upper left
  ah(end+1)=axes('Pos',[0.17 0.72 0.16 0.10],'tag','Spec01r');

  dd=fliplr(Idma.rho);
  if norm(dd(:,end))==0
     h=plot(dd(:,1:end-1),'ko-','Disp','w_n^{(D)}');
  else
     h=plot(dd,'o-');
     set(h(1),'Disp','w_n^{(D)}');
     set(h(2),'Disp','w_n^{(K)}');
  end
  sms(3); xtight; hold on; box on

  if exist('Inrg','var') && isfield(Inrg,'EE'), axis(axis); yl=ylim;
   % dd=0.2*Inrg.EE(1:min(end,128),:);
     dd=Inrg.EE;
        in=find(isnan(dd)); dd(in)=0;
        q=max(reshape(dd(:,end-2:end-1),[],1));
        i=find(sum(diff(dd,[],1).^2,2)>1E-2); dd(in)=nan;
     dd=dd(i,:);
     for j=1:2, i=j:2:size(dd,2);
        if j==1
             o={'Color',[.7 .7 .7]};
        else o={'Color',[1 .85 .7]}; end
        l=plot(i,dd(:,i)',o{:},'tag','eflow');
        mv2back(l);
     end
   % ylim([0 q/6]);
     ylim([0 2]);

     if isfield(Inrg,'NK') && min(size(Inrg.NK))==4
        n=int2str2(max(Inrg.NK)); o={'FontSize',10};
        postext(-0.05,1.40,'N_K^{}=%s (%s)',n{1:2},o);
        postext(-0.05,1.15,'N_D^{}=%s (%s)',n{3:4},o);
        if isfield(Inrg.param,'sym')
           s=Inrg.param.sym;
           s=regexprep(s,'SU(\d)','SU($1)');
           s=regexprep(s,',','\\otimes');
           s=regexprep(s,'(charge|spin|channel)*','_{$1}');
           postext(-0.05,1.70,'sym=%s',s,{o{:},'Color',[0 .5 0]});
        end
     end
  end

  xl=xlim;
  if xl(2)<60, set(gca,'XTick', [1 10:10:xl(2)]); end

% xlabel('NRG iteration');
  
  if exist('NKEEP','var')
     if exist('Inrg','var'), dd=Inrg.NK(:,1); else dd=NKEEP; end
     yl=ylim;
     plot(dd * (0.8*yl(2)/max(dd)),'r'); % mv2back(h(end));
     t=title(sprintf('N_{keep}^{max}=%g',max(dd)));
     set(t,'color',[.7 0 0]);
     xtight(1.05);
  end

  legdisp('Location','NorthEast','-detach','dy',-0.2,'xsc',0.7);

% --------------------------------------------------------------------- %
% --------------------------------------------------------------------- %

