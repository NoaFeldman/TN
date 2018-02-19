function set_global()
% Function: set_global()
% Wb,Aug07,07

% do not used variables as they are easily cleared
% e.g. using clr

  s=get(0,'UserData');

  if ~isstruct(s) || ~isfield(s,'err_count')
     s.err_count=struct('err',0,'wrn',0,'tst',0);
  end

  s.mroot=getenv('MYMATLAB');
  set(0,'UserData',s);

end

