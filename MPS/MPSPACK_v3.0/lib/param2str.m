function s=param2str(param,varargin)
% Function: param2str - convert parameter structure to string
% Usage: s=param2str(param [, OPTS, replace strings])
%
%    param is the parameter structure
%    replace strings are given by ..., 'str1', 'str1r', 'str2', 'str2r', ...
%    where str1 would be replace by str1r and so on.
%
% Options
%
%   '-tex'   escale some Greek letters and underscore
%   '-x',..  exclude fields matching grep pattern
%   'sep',.. field separator (', ')
%   'fmt',.. format for numbers and vectors
%
% Wb Jan 2007

  getopt ('init', varargin);
     tex = getopt('-tex'    );
     xpat= getopt('-x',   []);
     sep = getopt('sep',', ');
     fmt = getopt('fmt',  {});
  varargin=getopt('get_remaining');

  if isempty(param), s=''; return; end

  fn=fieldnames(param)'; m=length(fn);
  fn(4,:)={sep}; fn(2,:)={'='};

  if ~isempty(fmt), fmt={'fmt',fmt}; end

  for i=1:m
     if ~isempty(xpat) && ~isempty(regexp(fn{1,i},xpat))
     fn(:,i)={''}; continue; end

     d=getfield(param,fn{1,i});
     if ischar(d)
          fn{3,i}=[ '''' d '''' ];
     elseif ~isnumeric(d), fn(:,i)={''}; continue
     elseif numel(d)>4
          fn{3,i}=sprintf('[%dx%d array]',size(d));
     elseif length(d)>1
          fn{3,i}=[ '[' vec2str(d,fmt{:}), ']'];
     else fn{3,i}=vec2str(d,fmt{:}); end
  end

  for i=m:-1:1
     fn{4,i}='';
     if ~isempty(fn{1,i}), break; end
  end

  s=cat(2,fn{:});

% replace 
  if tex
     s=regexprep(s,'\\*\<(Lambda|Gamma|Delta|alpha|delta|omega)\>','\\$1');
     s=strrep(s,'_','\_');
     s=strrep(s,'epsd','\epsilon_d');
  end

  for i=2:2:length(varargin)
     s=strrep(s,varargin{i-1},varargin{i});
  end

end

