function banner(varargin)
% Function banner([flag,] str ...)
%
%     print text (handed over to sprintf) as banner
%     with flag setting the flashiness:
%       1   most flashy
%       2   intermediate
%       3   less flashy
%       4   least flashy
%      '%'  show as matlab comment with leading % (otherwise similar to 4)
%
% NB! info string is also written to window title

  if nargin==0
  eval(['help ' mfilename]); return; end

  persistent tlast

% flag=2-isbatch;

  if nargin>1 && isscalar(varargin{1}) %&& isnumeric(varargin{1})
       flag=varargin{1}; varargin=varargin(2:end);
  else flag=1; end

  nc=getcols(); if nc<74, nc=74; end

  switch flag
    case 1
    % also show elapsed time since last call to banner // Wb,Sep19,11
      if ~isempty(tlast)
           s=now-tlast; if s<1
                s=sprintf('  [%s]',datestr(s,   'HH:MM:SS'));
           else s=sprintf('  [%s]',datestr(s,'dd-HH:MM:SS')); end
      else s=''; end
      line=repmat('*',1,nc); tlast=now;
      fprintf(1,'\n\n%s\n\n  %s/%s\n  %s%s\n\n  ',...
      line, hostname, pwd, datestr(now),s);
    case 2
      line=repmat('*',1,nc); % tlast=[];
      fprintf(1,'\n\n%s\n\n  ',line);
    case 3
      line=repmat('=',1,nc); % tlast=[];
      fprintf(1,'\n%s\n  ',line);
    case 4
      line=repmat('-',1,nc); % tlast=[];
      fprintf(1,'\n%s\n  ',line);
    case '%'
      line=repmat('-',1,nc); % tlast=[];
      fprintf(1,'\n%% %s\n%% ',line);
    otherwise, % tlast=[];
  end

% eval(['! cd ~/Source/FigLet; figlet ' sprintf(varargin{:}) ]);
  if numel(varargin)>1
       s=sprintf(varargin{:});
  else s=cat(2,varargin{:}); % e.g. avoid interpretation of %
  end

  if flag~='%'
   % replace newlines with '\n  ' (indentation)
     s=regexprep(s,'\n([^#%])','\n  $1');
     fprintf(1,'%s\n',s);
  else
   % replace newlines with '\n  ' (indentation)
     i=' % '; i(1)=10;
     fprintf(1,'%s\n',strrep(s,i(1),i));
  end

  switch flag
    case 1, fprintf(1,'\n%s\n\n',line);
    case 2, fprintf(1,'\n%s\n\n',line);
    case 3, fprintf(1,'%s\n\n',line);
    case 4, fprintf(1,'%s\n\n',line);
    case '%', fprintf(1,'%% %s\n\n',line);
  end

  if flag==1, inl 1; end

% write info also to window title
% if ~isbatch && flag<=2
%  % NB! might fail if `s' contains carriage return or backslashes!
%    s=strrep(s,'\','');
%  % s=strrep(s,char(10),'; ');
%    i=find(s==char(10)); if ~isempty(i), s=s(1:i(1)-1); end
%
%    try 
%        eval(['! echo -ne "\033]0;matlab @  ${HOST}  ' ...
%        regexprep(repHome(pwd),'\$','\\$') ' : ' s '\007"'])
%    catch
%    end
% end

end

