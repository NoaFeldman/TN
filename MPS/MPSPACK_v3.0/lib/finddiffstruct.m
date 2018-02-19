function e=finddiffstruct(a,b,varargin)
% function e=finddiffstruct(A,B)
% 
%    Find fields that differ in between
%    structures or objects A and B.
% 
% Options
% 
%    '-v'  verbose mode
%    '-l'  recursive mode (enter fields that contain structs if different)
% 
% Wb,Aug08,16

  getopt('init',varargin);
     vflag=getopt('-v');
     rflag=getopt('-r'); % recursive in nested structures
     pre=getopt('--prefix','');
  getopt('check_error');

  if nargin<2 || nargout>1
     helpthis, if nargin || nargout
     error('Wb:ERR','invalid usage'), end, return
  end

  if isobject(a), a=struct(a); end
  if isobject(b), b=struct(b); end

  if ~isstruct(a) || ~isstruct(b)
     helpthis, if nargin || nargout
     error('Wb:ERR','invalid usage'), end, return
  end

  fa=fieldnames(a);
  fb=fieldnames(b); nl=repmat(' |  ',1,numel(find(pre=='.')));

  ff=intersect(fa,fb); e=0;

  if numel(ff)~=numel(fa), e=e+1;
     xa=setdiff(fa,fb); s=sprintf(', %s',xa{:});
     fprintf(1,'%s-> got %g extra field(s) in A: %s.%s\n',nl,pre,s(3:end));
  end
  if numel(ff)~=numel(fb)
     e=e+1; if e==1, fprintf(1,'\n'); end
     xb=setdiff(fb,fa); s=sprintf(', %s',xa{:});
     fprintf(1,'%s-> got %g extra field(s) in B: %s.%s\n',nl,pre,s(3:end));
  end

  if rflag
     o={'-r'}; if vflag, o{end+1}='-v'; end
  end

  for i=1:numel(ff)
     xa=getfield(a,ff{i});
     xb=getfield(b,ff{i});
     if ~isequal(xa,xb)
         e=e+1; if e==1, fprintf(1,'\n'); end
         if rflag && (isobject(xa) || isstruct(xa)) ...
                  && (isobject(xb) || isstruct(xb))
            if vflag, fprintf(1,['%s-> got different field %s.%s ' ...
               'that contains structure ...\n'],nl,pre,ff{i});
            end
            finddiffstruct(xa,xb,'--prefix',[pre '.' ff{i}],o{:});
         else
            fprintf(1,'%s-> got different field %s.%s\n',nl,pre,ff{i});
         end
     end
  end

  if vflag && ~e
     fprintf(1,'\n   got matching fields.\n');
  elseif e, fprintf(1,'\n');
  end

% wblog('TST',''); keyboard

end

