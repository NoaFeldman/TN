function label(varargin)
% Function: label() - adds x/y label and title to current figure
% Usage: label ([ah,] 'xstr', 'ystr', 'tstr' [, yorientation]);
%
%    If axes handles ah are present, apply labeles to every axis.
%
% Wb,Feb28,01

  if numel(varargin) && ~isempty(varargin{1}) && all(isaxis(varargin{1}))
       aflag=1; ah=varargin{1}; varargin=varargin(2:end);
  else aflag=0; ah=gca; end

  roty=0;
  if numel(varargin)
     if isnumeric(varargin{end})
        roty=varargin{1}; varargin=varargin(1:end-1);
     elseif isequal(lower(varargin{end}),'roty')
        roty=1; varargin=varargin(1:end-1);
     end
  end

  narg=numel(varargin);
  if narg<1 || roty && narg<2 || narg>3
     eval(['help ' mfilename]);
     if nargin || nargout, error('Wb:ERR','invalid usage'), end, return
  end

  tags = { 'XLabel','YLabel','Title' }; lh=cell(1,3);

  for i=1:narg
   % if input is [], then current/old label survives
     lh{i}=label_set(ah,tags{i},varargin{i});
  end

  if roty, set(lh{2},'Rotation',0); end

end

% -------------------------------------------------------------------- %
% -------------------------------------------------------------------- %
% if str is empty, then current/old label survives

function lh=label_set(ah,tag,str)

% isequal('',[]) yields 1 (!!)
  lh=get(ah,tag); if isempty(str) && ~ischar(str), return; end

% if label was not set yet, the fontsize is changed to default
% even though fontsize may have been set manually in the meantime !??
% fs=get(lh,'FontSize');

  if length(ah)>1, lh=cat(1,lh{:}); end
  set(lh,'String',str); % ,'FontSize',fs

end

% -------------------------------------------------------------------- %
% -------------------------------------------------------------------- %
