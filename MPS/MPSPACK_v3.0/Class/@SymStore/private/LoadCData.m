function q=LoadCData(sym,varargin)
% function q=LoadCData(sym,varargin)
%    
%    load specified CData structure.
%
% Wb,Jan30,15

  D0=get_dir(sym,'CStore');

  qstr=[varargin{:}];
  qstr=regexprep(qstr,'.*(','');
  qstr=regexprep(qstr,').*','');

  qdir=regexprep(qstr,'(\d+)\*,*','-');
  qdir=regexprep(qdir,'(\d+),*','+');

  qstr=regexprep(qstr,',(\d+)\*',';$1','once');
  qstr=regexprep(qstr,'\*','');

  f=[D0 '/' qdir '/(' qstr ').cgd'];
  if exist(f,'file') q=load2(f,'-mat');
  else error('Wb:ERR','\n   ERR file not found "%s"',f); end

end

