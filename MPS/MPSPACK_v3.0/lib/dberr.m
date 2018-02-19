function dberr(varargin)
% Wb,Jun29,11

  cmd='dbclear all ; dbstop if error';
  if nargin, cmd=[cmd, ' ; dbstop ' sprintf(' %s',varargin{:})]; end
% cmd

  evalin('caller',cmd);

end

