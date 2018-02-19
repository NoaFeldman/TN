function dispstack(S,k)
% Function dispstack(S [,k])
% Wb,Aug24,07

  if nargin<1 || nargin>2
     eval(['help ' mfilename]);
     if nargin || nargout, error('Wb:ERR','invalid usage'), end, return
  end

  if isempty(S)
     fprintf(1,'\n(empty dbstack)\n\n');
     return
  end

  if nargin<2, k=-1; end

  if isobject(S), S=obj2struct(S); end % Wb,Aug09,16

  if ~isfield(S,'file')
     if isfield(S,'stack')
        if isfield(S,'message'),    fprintf(1,'\n  Message: %s\n', S.message); end
        if isfield(S,'identifier'), fprintf(1,'       id: %s\n',S.identifier); end
        S=S.stack;
     else error('Wb:ERR','invalid usage'); end
  end

% indicator for current frame // Wb,Apr24,13
% (see bt.m)
% if ~isfield(S,'curr'), S(1).curr=0; end
% for i=1:numel(S)
%     if isempty(S(i).curr), S(i).curr=0; end
% end

  n=numel(S); vflag=0; inl 1

  if vflag
     fprintf(1,'\nDBSTACK\n\n');
     for i=n:-1:1
         s=sprintf('%s:%d', S(i).file, S(i).line);
         fprintf(1,'%3d  %-40s %s\n', n-i+1,s,S(i).name);
     end
  else
     for i=1:n % n:-1:1 % 1:n % 
         f=S(i).file; s={' ','',''};
         if ~isempty(find(f=='@',1))
          % got matlab class; use matlab convention:
          % @class/method (as this can be used with vi)
            s{2}=regexprep(f,'.*/(@[^/]*/).*','$1');
         else q=''; end

         f=regexprep(f,'^.*/','');

         if i==k, s{1}='*'; end
         f2=S(i).name;
            if isequal(f2,f(1:end-2)), f2=''; else f2=['(' f2 ')']; end
            s{3}=sprintf('%s:%d',f,S(i).line);
         fprintf(1,' %2d%s  %s%-30s %s\n',i,s{:},f2);
       % sprintf('%s:%d', regexprep(S(i).file,'^.*/',''), S(i).line), ...
     end
  end

  inl 1

end

