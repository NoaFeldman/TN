function setfig(fs)
% function setfig(fs)
% Wb,Apr04,13

% adapted from setax()

  if nargin~=1 || nargout
     helpthis, if nargin || nargout
     error('Wb:ERR','invalid usage'), end, return
  end

  if ~ischar(fs)
     if ~isnumeric(fs) || numel(fs)~=1 || ~isfig(fs)
        error('Wb:ERR','\n   ERR invalid figure handle');
     end
  else
     h=findall(groot,'type','figure','tag',fs);
     if numel(h)~=1
        if numel(h)>1, error('Wb:ERR',...
          '\n   ERR %g figures found with tag ''%s''',fs);
        else, error('Wb:ERR',...
          '\n   ERR no figure with tag ''%s'' yet',fs); 
        end
     end
     fs=h;
  end

  set(0,'CurrentFigure',fs);

end

