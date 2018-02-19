function m=checkOM(A,varargin)
% function m=checkOM(A [,opts])
%
%    Check for outer multiplicity in given QSpace
%
% Wb,Feb15,14

  getopt('init',varargin);
     vflag=getopt('-v');
  getopt('check_error');

  n=numel(A); m=zeros(size(A)); nout=0;
  for k=1:numel(A)
     [qq,I,D]=uniquerows(cat(2,A.Q{:}));
     [dd,id,d2]=uniquerows(D');
     im=find(dd>1); if ~isempty(im)
       for j=1:numel(im), l=im(j);
          if vflag
              getsub(A,[I{id{l}}])
          end
          if ~nargout, nout=nout+1; if nout==1, fprintf(1,'\n'); end
              fprintf(1,['   got multiplicity of %g (%g times) ' ... 
             '[%g/%g records]\n'], dd(l), d2(l), dd(l)*d2(l), size(A.Q{1},1));
          end
       end
       if vflag, fprintf(1,'\n'); end
       if numel(im)
          if ~vflag, fprintf(1,'\n'); end
          m(k)=d2(im)*dd(im)-numel(im); if ~nargout, fprintf(1,...
           '   got total outer multiplicity of %g/%g records.\n\n', ...
               m(k), size(A.Q{1},1));
          end
       end
     end
  end

  if all(m(:)==0) && ~nargout
     fprintf(1,'\n   got no outer multiplicity\n\n');
     clear m
  elseif ~nargout, clear m
  end

end

