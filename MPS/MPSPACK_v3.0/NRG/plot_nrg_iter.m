function ee=plot_nrg_iter(varargin)
% function ee=plot_nrg_iter(varargin)
% 
%    plot HK and HT data for a given iteration.
%
% Usage
%
%   (1) plot_nrg_iter(q);          struct with entries q.HK and q.HT 
%   (2) plot_nrg_iter(HK,HT);      HK,HT as specified
%   (3) plot_nrg_iter([NRG,]k);    load data from file
% 
% Wb,Sep01,12

  usage=[]; istr='';

  if nargin && isstruct(varargin{1}) && ...
     isfield(varargin{1},'HK') && isfield(varargin{1},'HT'), usage=1;
     varargin={ varargin{1}.HK, varargin{1}.HT, varargin{2:end} };
  elseif nargin==2 && ...
     isQSpace(varargin{1}) && isQSpace(varargin{2}), usage=2;
  elseif nargin && nargin<=2 && isnumber(varargin{end}), usage=3;
     if nargin==1
        varargin={ [getenv('LMA') '/NRG/NRG'], varargin{:} };
     end
     try 
        if ischar(varargin{1})
           f=sprintf('%s_%02g.mat',varargin{1},varargin{2});
           q=load(f,'HK','HT');
           varargin={ q.HK, q.HT, varargin{3:end} };
        else usage=[]; end
     catch
        wblog('ERR','invalid NRG data (%s)',f);
        return
     end
     istr=untex(repHome(f));
  end

  if isempty(usage) || numel(varargin)>2
     helpthis, if nargin || nargout
     error('Wb:ERR','invalid usage'), end, return
  end

ah=smaxis(1,1,'tag',mfilename); header('%M'); addt2fig Wb
setax(ah(1,1)); box on

  k=0; kk=[]; dy=[]; yn=nan; nq=numel(varargin);
  for i=1:nq
     dd=sort(cat(2,varargin{i}.data{:})); n=numel(dd);
     plot(k+1:k+n, dd,'o-','color',getcolor(i)); hold on
     if n && i<nq
        [h,hx]=ymark(dd(end),'k--','istr',sprintf('E_n=%.5g',dd(end)));
        set(hx,'VerticalAl','top');
        p=get(hx,'Pos'); p(2)=p(2)-0.2; set(hx,'Pos',p);
        yn=dd(end);
     end
     if n && i>1
        [h,hx]=ymark(dd(1),'k--','istr',sprintf('E_1=%.5g',dd(1)));
        set(hx,'VerticalAl','bottom');
        p=get(hx,'Pos'); p(2)=p(2)+0.2; set(hx,'Pos',p);
        dy(i-1)=dd(1)-yn;
     end
     kk(end+1)=k+0.5; k=k+n;
     if nargout, ee{i}=dd; end
  end

  h=findall(gca,'type','text','tag','ymark'); n=numel(h);
  x=xlim; dx=diff(x)/n; x=dx/10:dx:x(2);
  for i=1:n
      p=get(h(i),'Pos'); p(1)=x(n-i+1); p(3)=10; set(h(i),'Pos',p);
  end
  set(h,'FontW','normal','FontS',10,'Margin',1E-3,'BackgroundC','w');
% keyboard

  sms(4); xtight;
  if numel(kk)>1
     title2('Nkeep=[%s]',vec2str(kk(2:end)-0.5));
     xmark(kk(2:end),'k:');
  end

  if ~isempty(istr), title(istr); end

end

