function v=mlver()
% function v=mlver()
% 
%    returns verion number based on release date
%
% Wb,Aug03,16

  v=lower(regexprep(version,'.*\(R20(.*)\).*','$1'));

  if isempty(regexp(v,'^\d+[a-z]')), version, v
     error('Wb:ERR','\n   ERR unexpected matlab version !?');
  end

  v=str2num(v(1:end-1)) + (v(end)-'a'+1)/10;

end

