function h=title2(varargin)
% function h=title2(varargin)
% 
%    like title() but also accepts sprintf syntax as in 
%    title(sprintf(...))
% 
% Wb,Feb06,10

  if nargin && isaxis(varargin{1})
     ah0=gca;
     setax(varargin{1}); varargin=varargin(2:end);
  else ah0=[]; end

  xflag=0;
  while numel(varargin) && ischar(varargin{1})
     switch varargin{1}
       case '-append',  xflag=+1;
       case '-prepend', xflag=-1;
       otherwise break;
     end
     varargin=varargin(2:end);
  end

  narg=numel(varargin);
  if ~narg
     if nargout, h=get(gca,'Title'); else, title(''); end
     return
  end

  if narg==1
     s=regexprep(varargin{1},'\\+','\\');
  else
     s=regexprep(varargin{1},'\\n([^a-z])',[char(10) '$1']);
     s=regexprep(s,'\\+','\\\\');
     s=sprintf(s,varargin{2:end});
  end

  if xflag
     s={ get(get(gca,'Title'),'string'), '; ', s};
     if xflag<0, s=s(end:-1:1); end
     s=cat(2,s{:});
  end

  title(s);

  if ~isempty(ah0), setax(ah0); end

end

