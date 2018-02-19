function [oo,aa,ah,I] = rsmoothSpec(Om,Aa,varargin)
% Usage: [oo,aa,ah,I] = rsmoothSpec(Om,Aa [,opts])
%
%   subsequent smoothening of spectral function
%
% Options:
%
%   Idma    if arg[3] is structure, T will be read out of it (if exists)
%   'T',..  temperature (for plot only)
%   'adisp',{..}  set of labels used as Disp in aa data
%   'afac'  factor when plotting smooth spectral function (pi*Gamma)
%   'nfac'  extra factor (1/nfac) when using SU(N) symmetries
%   'skip'  skip smooth data for |om|<skip (0)
%
%   'nofig' do not show data in figure
%   'xli'   limits for x axis on inset
%   'yli'   y-limits for main panels
%   'raw'   input data Aa is considered raw data
%   'RAW'   same as 'raw', but also plots raw data
%
% NB! remaining options are handed over to getSmoothSpec(Om,Aa,....)
%
% Wb,May12,06

  global N param

  if nargin<2, eval(['help ' mfilename]); return; end

  Idma=[];
  if length(varargin)>0
     if isstruct(varargin{1})
     Idma=varargin{1}; varargin=varargin(2:end); end
  end

  if ~isempty(param)
     TK=TKondo; xli=5*TK;
     if isfield(param,'B') && param.B>TK, xli=5*param.B; end
     if isequal(xli,0), xli=[]; end
  else
     xli=2E-4;
  end

  getopt ('INIT', varargin);
     T     = getopt('T',     []);
     eps   = getopt('eps',   []);
     reps  = getopt('reps',  []); % same as eps, put relative to T
     nofig = getopt('nofig'    );
     afac  = getopt('afac',  []);
     adisp = getopt('adisp', {});
     nfac  = getopt('nfac',   1);
     xli   = getopt('xli',  xli);  % xlim for inset
     yl    = getopt('yli',   []);  % ylim for main panels
     skip  = getopt('skip',   0);
     sigma = getopt('sigma', -1);
     keepA0= getopt('keepA0', 0);

     if getopt('raw'), rflag=1; 
     elseif getopt('RAW'); rflag=2; else rflag=0; end

  varargin=getopt('get_remaining');

  if length(xli)==1, xli=[-1 1]*xli; end

  if isempty(T) && isfield(Idma,'T') && ~isnan(Idma.T)
  T=Idma.T; end

  if sigma<=0 && ~isempty(param) && isfield(param,'Lambda')
     if param.Lambda==2, sigma=0.6; else
      % Lambda=2 => sigma=0.6
        sigma = log(param.Lambda) * (0.6/log(2));
     end
  end

  if isempty(eps) && ~isempty(T) && T>0
     if ~isempty(reps), eps=T*reps; else eps=T; end
  end
  if ~isempty(eps), varargin(end+1:end+2)={'eps',eps}; end

% i=length(Om)/2; a0=Aa(i:i+1,:); aref=1E3*norm(Aa(i+[-5:-2, 3:6],:));
  if keepA0, o={'-a0'}; else o={}; end
  if sigma>0, setopts(o,sigma); end

  if nfac~=1, Aa=Aa*(1/nfac); end

% wblog('TST','finished (0)'); pause(0);
% ************************************************** %
  [oo,aa,I] = getSmoothSpec(Om,Aa,varargin{:},o{:});
% ************************************************** %
% wblog('TST','finished (1)'); pause(0);

  if skip>0
     i=find(abs(oo)<skip);
     oo(i)=[]; aa(i,:)=[];
     I.skip=skip;
  end

  if nofig, ah=I; return; end

  if isempty(afac)
     if ~isempty(param)
        if isfield(param,'fJ') || isfield(param,'J')
           afac=(pi^2)/2; facstr='\pi^2/2';
        elseif isfield(param,'Gamma')
           afac=pi*max(abs(param.Gamma(:))); facstr='\pi\Gamma';
        else afac=1; facstr=''; end
     else
        wblog('WRN','param not available - scale max_A to 1.');
        afac=1/max(abs(aa(:))); facstr=sprintf('%.4g',afac);
     end
  else
     if iscell(afac)
        if numel(afac)~=2 || ~isnumeric(afac{1}) || ~ischar(afac{2})
        afac, error('Wb:ERR','invalid afac'); end

        facstr=afac{2}; afac=afac{1};
     elseif afac~=1
        facstr=sprintf('%.4g',afac);
     else facstr=''; end
  end

% afac=afac;

  I.afac={afac,facstr};
  I.nfac=nfac;

  if rflag
     a=sum(Aa); a=a(find(abs(a)>0.5 & abs(a-round(a))<1E-3));
     if ~isempty(a), a=a(1); end
  end

% wblog('TST','finished ...'); pause(0); 
% Matlab online session died here:
% 1 times // Wb,Jun20,13
% 1 times // Wb,Jun21,13
ah=smaxis(2,1,'tag',mfilename,'DY',0.02);
% wblog('TST','finished ...'); pause(0);
ah(3)=axes('Pos',[0.73 0.72 0.16 0.10]); % upper right inset
% wblog('TST','finished ...'); pause(0);

  if ~isempty(param), TK=TKondo; else TK=nan; end

  if isempty(yl)
     yl=[-0.05,  1.2]*afac*max(aa(:));
     if diff(yl)<=0, yl=[]; end
  end

% afac=1;

  lfmt0 = {'color', [.7 .7 .7],'tag','raw'};
  lfmtp = {'color', [.7 .7  1],'tag','raw'};
  lfmtn = {'color', [ 1 .7 .7],'tag','raw'};
  if rflag, lg1='cc raw data / \omega'; else lg1='A0'; end
  if rflag, lg2='cc raw data / \omega'; else lg2='A0'; end

  astr='A(\omega)';
  if ~iscell(adisp) && ischar(adisp)
     switch adisp
        case {'spin_ud','spin'}
        adisp={'A_{\uparrow}(\omega)','A_{\downarrow}(\omega)'};
        case 'spin_du'
        adisp={'A_{\downarrow}(\omega)','A_{\uparrow}(\omega)'};
     otherwise; end
  end
  if isempty(adisp), adisp={astr}; elseif ~iscell(adisp)
  adisp, error('Wb:ERR','\n   ERR invalid adisp'); end

  ox=Om(:); ox(find(ox==0))=nan; ip=find(ox>0); in=find(ox<0);

  oa={'tag','Adata'};

% wblog('TST','finished ...'); pause(0);
setax(ah(1))

  h2=plot(oo,afac*aa,oa{:}); hold on; if ~isempty(yl), ylim(yl), end
  fliphh(h2);

     n=min(numel(adisp),numel(h2));
     for i=1:n, set(h2(i),'Disp',adisp{i}); end

  if rflag>1
     ar=afac*Aa; or=Om;
   % for i=1:size(ar,2), ar(:,i)=ar(:,i) ./ abs(ox);  end
   % ar*(2/max(abs(a1(:))));
     ar=ar./repmat(abs(ox),1,size(Aa,2));

     nl=log(10)/abs(mean(diff(log(abs(ox(2:5))))));
     m=round(nl/250); % down sample to nlog of about 100
     if m>1
        l=mod(size(ar,1),m);
        if l, l=m-l;
           ar(end+(1:l),:)=0;
           or(end+(1:l))=or(end);
        end
        s=size(ar); ar=permute(sum(reshape(ar,[m,s(1)/m,s(2)]),1),[2 3 1]);
        s=length(or)/m; or=permute(sum(reshape(or,[m,s]),1),[2 1]);
     end

   % h1=plot(or,ar/max(ar(:)), lfmt0{:});
     h1=plot(or,ar,lfmt0{:});

     if ~isempty(I.a0), y=ylim;
     plot([0 0], [0 y(2)],'-', lfmt0{:},'LineW',2); end

     mv2front(h2);

     for i=1:length(h1)
        set(h1(i),'Color',0.7+0.2*getcolor(i));
     end
     set(h1(1),'Disp',lg1);
  else 
   % h1=plot(or,afac*ar, lfmt0{:});
  end

% keyboard 

  if ~isempty(afac) && isempty(a==1), a=[a 1]; end
% keyboard

  xlim([-1 1]*0.2);
  ymark([0 a],'k--');

  if ~isempty(facstr), astr=[astr '  \ast ' facstr]; end

  label('\omega',astr);
  legdisp -flip

% wblog('TST','finished ...'); pause(0);
setax(ah(3))

  h=plot(oo,afac*aa,oa{:});

    n=min(numel(adisp),numel(h));
    for i=1:n, set(h(i),'Disp',adisp{i}); end

  mv2front(h(1)); hold on
  
  if diff(xli)>0, xlim(xli); else
     if ~isempty(xli) && ~all(isnan(xli))
     wblog('WRN','invalid xli=[%s]',vec2str(xli)); end
     xtight
  end

  ytight('view'); l=ylim; ytight(1.2,'view'); l2=ylim; l2(1)=l(1); ylim(l2);
  ymark([0 a],'k--');

  if ~isempty(param) && isfield(param,'B') && param.B~=0
  xmark([-1 1]*param.B,'k--'); end

% keep default tag='xmark' => allows to use xmark -fix (!) %% Wb,Apr30,12
  ok={'Color',[.5 .75 .5],'LineW',2}; % ,'tag','TK'

  if isset('TK'), xmark([-1 1]*TK,ok{:}); end

% wblog('TST','finished ...'); pause(0);
setax(ah(2))

  if rflag>1
     jp=find(or>0); jn=find(or<0);
     h2=semilogx(abs(or(jn)),ar(jn,:), lfmtn{:}); hold on
     h1=semilogx(abs(or(jp)),ar(jp,:), lfmtp{:});
     set(h1(1),'Disp',lg1);
  else
   % h1=semilogx(abs(or),afac*ar);
  end

  ip=find(oo>0);
  in=find(oo<0);
  hn=semilogx(-oo(in), afac*aa(in,:),oa{:}); fliphh(hn); hold on
  h2=semilogx( oo(ip), afac*aa(ip,:),oa{:}); fliphh(h2);

  set(hn,'LineSt','--'); set(hn(1),'Disp','(\omega<0 data)');
% set(hn,'tag','omega<0'); tag already set to Adata!

    n=min(numel(adisp),numel(h2));
    for i=1:n, set(h2(i),'Disp',adisp{i}); end

  mv2front(h2); hold on; xtight

% if ~isempty(yl), set(gca,'ylim',yl); end
  if ~isempty(yl), ylim(yl); else
  ylim(get(ah(1),'YLim')); end

  topt = { 'HorizontalAlignment','Center', ...
  'VerticalAlignment','middle' };

  if ~isempty(param) && isfield(param,'Lambda') && ~isempty(N)
     xl=param.Lambda^(-N/2);
     lh=plot([xl xl], ylim, 'LineWidth',4, 'color', [1 .95 .5]);
     text(xl,.8, sprintf('\\Lambda=%g, N=%g', param.Lambda, N), topt{:});
     mv2back(lh);
  end

% yp=ylim; yp=yp(1)+[0.4 0.2]*diff(yp);

  if ~isempty(param)
     xmark(TK,'top','istr','T_K',ok{:});
  end

  if isfield(I,'eps')
     c2=[1 1 1]*0.8;
     xmark(I.eps,'top','istr','\epsilon','Color',c2,'LineW',3);
  end

  if ~isempty(T)
  xmark(T,'top','istr','T'); end
  ymark([0 a],'k:');

  s={};
  if isfield(I,'sigma')
     s{end+1}=sprintf('\\sigma=%.3g',I.sigma);
   % text(.85,.75,sprintf('\\sigma=%.3g',I.sigma),'Units','norm',...
   %  'EdgeColor','k','BackGroundColor','w','margin',10,'tag','sigma');
  end

  if isfield(param,'ALambda') && ~isempty(param.ALambda)
     s{end+1}='using A_\Lambda';
  end
  if ~isempty(s)
     s=sprintf('(%s)',strhcat(s{:},'sep','; '));
     text(.88,.70-0.05*numel(adisp),s,'Units','norm','HorizontalAl','center');
  end

  if ~isempty(I.a0)
     for i=1:size(I.a0,2)
         s=sum(I.a0(:,i)); if sum(abs(s))<1E-10, continue; end
         postext(.04,1-i*.1, ...
         'A_%d: %.6g \\delta(\\omega)',i,s, {'Color',get(h2(i),'Color')});
     end
  end

  label('|\omega|',astr);
  legdisp -flip

  set(ah(1),'tag','Spec01');
  set(ah(2),'tag','Spec02');
  set(ah(3),'tag','Spec01i');

% set(gcf, 'CurrentAxes', ah(1));
  if nargout==0, clear oo aa; end

% wblog('TST','finished (2)'); pause(0);

end

