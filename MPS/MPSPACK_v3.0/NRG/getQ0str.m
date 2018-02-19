function Qs=getQ0str(H,varargin)
% function Qs=getQ0str(H)
% Wb,May30,11

  getopt('init',varargin);
     istr=getopt('istr','Q');
     eps =getopt('eps',0.01);
  getopt('check_error');

  if isfield(H,'HK') && isfield(H,'E0')
     Inrg=H;
     l=numel(Inrg.E0)-2;
     e0=Inrg.E0(l:l+1); if e0(1)>e0(2), l=l+1; end
     H=Inrg.HK(l);
  end

  if isdiag(QSpace(H))<2
     [ee,H]=eig(QSpace(H)); % Wb,Apr26,13
  end

  dd=H.data; e0=min(cat(2,dd{:}))+eps; n=numel(dd); Qs={}; im=[];
  for i=1:n
     if ~isempty(find(dd{i}<=e0,1))
        im(end+1)=i; s=sprintf(' %g',H.Q{1}(i,:));
        Qs{end+1}=sprintf('[%s]',s(2:end));
     end
  end

  if numel(Qs)>1
       Qs=sprintf(',%s',Qs{:}); Qs=sprintf('\\{%s\\}',Qs(2:end));
  else Qs=Qs{1}; end

  if ~isempty(istr), Qs=sprintf('%s=%s',istr,Qs); end

end

