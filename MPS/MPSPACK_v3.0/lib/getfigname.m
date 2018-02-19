function fname=getfigname(varargin)
% function fname=getfigname([fh])
% Wb,Apr05,11

  f=[];

  while numel(varargin), v=varargin{1};
     if ishandle(v) && isequal(get(v,'type'),'figure')
        f=v; varargin=varargin(2:end);
     else
        helpthis, if nargin || nargout
        error('Wb:ERR','invalid usage'), end, return
     end
  end

  if isempty(f), f=gcf; end

  fname=getuser(f,'fname');
  if isempty(fname)
     fname=getuser(f,'mat');
     if isempty(fname)
        fname=get(f,'tag');
        if isempty(fname), fname='figure'; end
     else
        [p,n,x]=fileparts2(fname); t='';
      % t=get(f,'tag'); if ~isempty(t), t=[t '_']; end
        fname=[p t n];
     end
  end

end

