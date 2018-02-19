function setglobal(varargin)
% Function setglobal('var1','var2',...)
%
%   set input arguments global in current workspace
%
% Wb,Feb02,05

% NB! isglobal is an obsolete function and will be
% discontinued (see help isglobal) use this instead:
% if ~isempty(whos('global','variable')), ... 

% NB! iff local and global instance coexist
% clear global does not clear local instance

% fflag = force local instance if it exists global without warning
  fflag=0; if isequal(varargin{1},'-f')
  fflag=1; varargin=varargin(2:end); end

% use sgFlags to communicate between this function and calling function
  if ~isempty(whos('global','sgFlags'))
  wblog('WRN','sgFlags already exists'); end;
  
  evalin('caller','global sgFlags'); global sgFlags

  nv=length(varargin);
  for iv=1:nv, vn=varargin{iv};
  % il: exists local (and maybe also global if il.global)
  % ig: exists global
  
    evalin('caller', sprintf(...
    'sgFlags.il=whos(''%s''); sgFlags.ig=whos(''global'',''%s'');', vn,vn));

  % easy stuff first: variable not even local yet or already global
    if isempty(sgFlags.il)
    evalin('caller',['global ' vn]); continue;
    elseif sgFlags.il.global, continue; end

  % from here on, *local* instance exists
    if ~isempty(sgFlags.ig) && ~fflag
     % if variable coexists both local and global
     % *** ALWAYS PICK LOCAL ***
       fprintf(1,['\nWarning: %s - variable %s coexists ' ...
       '[local] and global!!\nWarning: %s\n\n'],mfilename,vn,lineno('all'));
  % else no global instance yet
    end

  % avoid warning message by creating a local copy of <variable>
    evalin('caller', sprintf(...
    'var__=%s; clear %s; global %s; %s=var__; clear var__',...
    vn,vn,vn,vn));

  end

  clear global sgFlags

end

