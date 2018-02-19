function Iw=getrhoESpectra(R,varargin)
% function Iw=getrhoESpectra(RR,varargin)
%
%    graphical analysis of entanglement spectrum of density matrix,
%    splitting the data into symmetry subspaces.
%
% alternative usage: Iw=getrhoESpectra(IR,kk,varargin)
%
%    where IR specifies a larger set, where only elements kk are
%    to be considered. IR is either of the following ...
%     1) RR        set of density matrices
%     2) IR[.R]    use output structure from getrhoNRG_resall.m
%     3) <nrgtag>  (string that locates NRG data) => use EK data instead of Rho
%     4) { (RR|IR, <nrgtag> } show both, energy and density matrix spectra
%
% see also getrhoNRG_anal
% Wb,Apr18,11

  if nargin>1 && isnumeric(varargin{1}) && isvector(varargin{1})
  gotkk=1; else gotkk=0; end

  if iscell(R)
   % usage 2
     if ~gotkk, error('Wb:ERR','\n   ERR invalid usage'); end
     mixed_usage(R,varargin{:});
     return
  end

  if gotkk
   % usage 2 => got kk set
     kk=varargin{1}; varargin=varargin(2:end);
     [R,I]=select_data(R,kk); nR=numel(R); 
   else 
   % usage 1
     [R,I]=select_data(R); nR=numel(R); 
   end

   getopt('init',varargin);
      yl = getopt('yl',[]); % ylim;
      y2 = getopt('y2',[]); % range for show_deg
      IST= getopt('istr',I.IST);
      eps= getopt('eps',0.1);
   getopt('check_error');

   odeg={ 'eps',eps,'dx',[0.18 0.35],'n',256,'fs',9 };
   % 'BackgroundC','w','Margin',1E-3

   o={}; switch nR
      case 1, o={'fpos',[890 635 710 485]};
      case 2, o={'fpos',[890 390 710 730],'dy',0.02};
   otherwise, o={'fpos',[890 70 710 1050]}; end

if ~I.Eflag, t=mfilename; else t=[mfilename '_E']; end
ah=smaxis(nR,2,'tag',t,'dx',0.02,o{:},'Y0',0.02);
evalin('caller','nrg_header'); addt2fig Wb
set(gcf,'defaultAxesLineWidth',1.5);

   if gotkk
        header('%M :: k=[%s]%s',vec2str(kk),I.istr);
   else header('%M :: %s',I.istr); end

   for i=1:nR
      if i<=numel(IST), istr=IST{i}; else istr=sprintf('k=%g',kk(i)); end
      getrhoESpectra_aux1(ah(i,:),R(i),I.Eflag,yl,y2,istr,odeg);
   end

   for h=reshape(ah(1:end-1,:),1,[])
      set(h,'XTickLabel',[]);
      set(get(h,'xlabel'),'string',''); 
   end

end

% -------------------------------------------------------------------- %
% -------------------------------------------------------------------- %

function mixed_usage(RC,kk,varargin)

   if numel(RC)~=2, error('Wb:ERR','\n   ERR invalid usage (cell R)'); end

   getopt('init',varargin);
      yl = getopt('yl',[]);
      y2 = getopt('y2',[]);
      eps= getopt('eps',0.1);
   getopt('check_error');

   global param

   if numel(eps)==1, eps=[eps/5 eps]; end
   odeg={'dx',[0.18 0.35],'n',256,'fs',9};
   odegE={'eps',eps(1), odeg{:}};
   odegR={'eps',eps(2), odeg{:}};

   RC=reshape(RC,[],1); nRC=numel(RC); nR=numel(kk); 
   for j=1:nRC
      [R{j},I(j)]=select_data(RC{j},kk);
   end
   kk=kk-1;

   e=cat(1,I.Eflag); t=cell(size(e));
   t(find( e))={'E'};
   t(find(~e))={'R'}; t=cat(2,t{:}); % ER ER2

   for ik=1:nR
      if nR>1, tag=sprintf('%s%g',t,ik); else tag=t; end
      mat=sprintf('%s_%s_D%g',mfilename,tag,param.D);
      if isfield(param,'U')
          if norm(param.U+2*param.epsd)>1E-12
          mat=sprintf('%s_ed%02g',-10*param.U/param.epsd); end
      end
      if isfield(param,'B') && param.B
          mat=sprintf('%s_B%02g',param.B);
      end

      tag=[mfilename,'_',tag];
      
      ah=smaxis(2,nRC,'tag',tag,'dx',0.02,'X0',-0.02,'Y0',0.02,...
        'fpos',[970 480 630 630],'dy',0.08,'zoom',0.85);
      set(gcf,'defaultAxesLineWidth',1.5); mvaxis(ah(1,:),[0 0.02]);
      evalin('caller','nrg_header'); addt2fig Wb
      header('%M :: k=%g \t [%s]; %s',kk(ik),vec2str(kk),I(1).istr);

      setuser(gcf,'mat',mat);

      for j=1:nRC
         if I(j).Eflag, odeg=odegE; else odeg=odegR; end
         getrhoESpectra_aux1(ah([2 1],j),R{j}(ik),I(j).Eflag,yl,y2,'',odeg);
      end

    % adjust labels, ticks, ...
      for h=reshape(ah(:,2:end-1),1,[])
         set(h,'YTicklabel',[]); ylabel(h,'');
      end
      for j=[1 nRC]
         set(ah(2,j),'YTickLabelMode','auto');
         ylabel(ah(1,j),get(get(ah(2,j),'YLabel'),'String'));
      end
      for h=ah(:,end)'
         set(h,'YAxisLoc','right');
         set(get(h,'YLabel'),'Rotation',-90,'VerticalAl','bottom');
      end

      ht=ah';
      for i=1:numel(ah)
         setax(ht(i)); poslabel(i,'NW','dx',[0 -0.02],'FontS',15);
      end

    % drawnow
    % if ik<nR, if ik>1, pause, else wbpause; end, end
   end

end

% -------------------------------------------------------------------- %
% -------------------------------------------------------------------- %

function [R,I]=select_data(R,kk)

   if nargin<2, kk=1:numel(R); end
   nR=numel(kk); 

   if ~ischar(R), R=R(kk); end

 % NB! RR(k).R => RR(k).k1 = k1-1 => offset by 1 (!)
 % eg. assuming output from getrhoNRG_resall
   kk=kk-1;

   I=struct('Eflag',0,'istr','','IST',{{}});

   if ischar(R), I.Eflag=1; nrg=R; R=[];
   % use energy data instead of input density matrices: R = NRGDATA
   % for comparison with entanglement spectra
     for k=kk
        q=load(sprintf('%s_%02g.mat',nrg,k-1));
        q=diag(QSpace(q.HK));
        if ~isempty(R), R(end+1)=q; else R=q; end
     end
   end

   if nargin>1
      getbase Inrg
      if ~isempty(Inrg), s=vec2str(Inrg.E0(kk),'fmt','%.4g');
           I.istr=sprintf(' @ {\\Delta}E_0=[%s]',s);
      else I.istr='(Inrg not available)'; end
      s=cell(size(kk));
      s(find(mod(kk,2)==0))={'even'};
      s(find(mod(kk,2)==1))={'odd'}; I.IST=s;
   end

   if isstruct(R) && isfield(R,'R')
    % eg. assuming output from getrhoNRG_resall
      R=cat(2,R.R);
   end
end

% -------------------------------------------------------------------- %
% -------------------------------------------------------------------- %

function EE=getrhoESpectra_aux1(ah,R,Eflag,yl,y2,istr,odeg)

   Q=R.Q; if numel(Q)~=2 || ~isequal(Q{:})
      error('Wb:ERR','\n   ERR invalid density matrix R'); end

   [rr,I]=eigQS(R); rr=rr(:,1);
      EK=I.EK; e=norm(rr(find(rr<=0))); if e>1E-12
      error('Wb:ERR','\n   ERR invalid density matrix (%.3g)',e); end

   Q=round(1E3*I.EK.Q{1})*1E-3; [N,n]=size(Q); II=zeros(N,n);
   q=cell(1,n); I=q; D=q; t='';
 % N: total number of symmetry sectors
 % n: total number of symmetry labels
   for i=1:n
     [q{i},I{i},D{i}]=uniquerows(Q(:,i)); m=numel(I{i});
     ii=zeros(N,1); for j=1:m, ii(I{i}{j})=j; end
     II(:,i)=ii;
   end

   if n<2, II(:,end+1:2)=1; end

 % EE: eigenspectrum resolved w.r.t. first 2(!) symmetries only
 % ==> EE(i1,i2,:) energies for symmetry q{1}(i1) and q{2}(i2)
   for i=1:numel(EK.data)
     dd=EK.data{i}; 
     if Eflag==0 % got rho => entanglement spectrum
        dd(find(dd<=0))=[]; dd=sort(-log(dd));
     end
     m=numel(dd); dd(find(dd==0))=-1; % tags: MARK_ZERO
     EE(II(i,1), II(i,2), 1:m)=dd;
   end

   EE(find(EE==0))=nan;
   EE(find(EE==-1))=0; % tags: MARK_ZERO

   if 0
      n=size(EE,3); setax(ah(1)); 
      k=1:n; while keyiter(k)
         mp(EE(:,:,k),'cmap','-gca'); header('k=%g/%g',k,n);
      end
   end

   if size(Q,2)==2
      qlabel={'charge Q', 'spin 2S_z'};
   else % tags: Q1 Q2 Q_1 Q_2
      for i=1:size(Q,2), qlabel{i}=sprintf('Q_{%g}',i); end
   end

   r2=0.40; lo={'Color',[.7 .7 .7]}; s=size(EE);

 % a=(5/6);
   a=1;
   if isempty(yl) % ylim
      if Eflag, yl=[-0.5 6]*a;
      else      yl=[-2.2  27]*a; end
   end

   if isempty(y2) % range for show_deg
      if Eflag, y2=4*a;
      else      y2=21*a; end
   end

   if Eflag, dy=0.5;
   else      dy=1.5; end

   odeg={odeg{:},'y2',y2,'dy',dy};

setax(ah(1)); 

   if 1 && ~Eflag
    % shifted entenglement spectrum
      EE=EE-min(EE(:)); t=' (shifted)';
   end

   dd=reshape(EE,s(1),[]); % dd(find(dd==0))=nan;
 % plot(q{1},dd,'+','MarkerS',8);

   for i=1:size(dd,1), x=q{1}(i);
    % avoid accessive number of handles (slow!) => break lines by nan's
      xx=repmat([x-r2;x+r2],1,size(dd,2)); xx(3,:)=nan;
      yy=repmat(dd(i,:),3,1);
      plot(xx(:),yy(:),lo{:}); hold on
      show_deg(dd(i,:),gca,'xt',x,odeg{:}); % postext
   end

   xl=xtight(1.02); xt=ceil(xl(1)):xl(2); set(gca,'XTick',xt);
 % yl=ytight(1.02,'y1',-1);
   ylim(yl); if yl(1)<0, ymark(0,'k:'); end

   xlabel(qlabel{1})
   if ~isempty(istr), istr=['   (' istr ')']; end
   if Eflag
        ylabel(['rescaled energy' istr]);
   else ylabel(['entanglement spectrum' istr t]); end
 % title(istr);

setax(ah(2))

   dd=reshape(permute(EE,[2 1 3]),s(2),[]); % dd(find(dd==0))=nan;
   plot(q{2},dd,'+','MarkerS',8);

   for i=1:size(dd,1), x=q{2}(i);
    % avoid accessive number of handles (slow!) => break lines by nan's
      xx=repmat([x-r2;x+r2],1,size(dd,2)); xx(3,:)=nan;
      yy=repmat(dd(i,:),3,1);
      plot(xx(:),yy(:),lo{:}); hold on
      show_deg(dd(i,:),gca,'xt',x,odeg{:}); % postext
   end
   xlim(xl); ylim(yl); ymark(0,'k:');
   set(gca,'XTick',xt); % ,'YTickLabel',[]

   xlabel(qlabel{2}); %,'entanglement spectrum \xi_i');

 % keyboard

end

% -------------------------------------------------------------------- %
% -------------------------------------------------------------------- %

