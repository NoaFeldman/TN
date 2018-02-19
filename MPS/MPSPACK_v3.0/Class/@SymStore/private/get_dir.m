function D=get_dir(sym,varargin)
% function D=get_dir(sym,varargin)
% Wb,Jan30,15

% sym=regexprep(sym,'(S[pU])\(*(\d+)\)','$1$2\n');
  D=[getenv('RC_STORE') '/' sym];

  if ~exist(D,'dir'), error('Wb:ERR',...
     '\n   ERR RCStore/%s does not yet exist',sym); end

if nargin<=1, return; end

  for i=1:numel(varargin)
     if ~ischar(varargin{i}), error('Wb:ERR','\n   ERR invalid usage'); end
  end

  D=[D sprintf('/%s',varargin{:})];

  if ~exist(D,'dir'), error('Wb:ERR',...
    '\n   ERR got non-existent directory\n   ERR %s',D); end

end

