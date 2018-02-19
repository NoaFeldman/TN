function dd=cat2(ic,varargin)
% function dd=cat2(d,varargin [,{value}])
% similar to cat, but adjust dimensions if required
%
%    last argument if of cell type {value} sets value
%    for space added; default: 0
%
% Wb,Jul23,09

  if ~nargin || ~any(ic==(1:3))
     eval(['help ' mfilename]);
     if nargin || nargout, error('Wb:ERR','invalid usage'), end, return
  end

  x=0; % default padding value
  if iscell(varargin{end}) && numel(varargin{end})==1 && ...
     isnumeric(varargin{end}{1})
     x=varargin{end}{1}; varargin=varargin(1:end-1);
  end

  n=length(varargin); if n<1, dd=[]; return; end

  s=cell(1,n); rr=zeros(1,n);
  for i=1:n
     s{i}=size(varargin{i}); rr(i)=length(s{i});
  end

  r=max(ic,max(rr));
  ss=ones(n,r); % adding singletons if required
  for i=1:n, ss(i,1:rr(i))=s{i}; end

  k=1:r; k(ic)=[];

  smin=min(ss(:,k),[],1);
  smax=max(ss(:,k),[],1);

  if ic<=2 && r>2 % code below assumes matrizes!
  error('Wb:ERR','invalid usage'); end

  if smin~=smax
     switch ic
       case 1, for i=1:n, varargin{i}(:,(end+1):smax)=x; end
       case 2, for i=1:n, varargin{i}((end+1):smax,:)=x; end
       case 3, for i=1:n, varargin{i}((end+1):smax(1),(end+1):smax(2),:)=x; end
       otherwise, error('Wb:ERR','invalid usage (ic=%g)',ic);
     end
  end

  dd=cat(ic,varargin{:});

end

