function info(A,s,cflag)
% function info(A [,istr,cflag])
%
%    header routine used in display()
%    Options: istr (info string) and cflag (compact info).
%
% Wb,Mar01,08

  if nargin<3, cflag=0; end
  if nargin<2, nm=inputname(1);
     if ~isempty(nm), s=sprintf('%s = ',inputname(1));
     else s=''; end
  elseif ~ischar(s)
     helpthis, error('Wb:ERR','invalid usage (istr must be string)')
  end

  nA=numel(A); if nA>1
     if isempty(s), nm='ans'; else nm=regexprep(s,'[= ]*$',''); end
  end

  for k=1:nA % Wb,Jun04,16
     if nA>1, s=sprintf('%s(%g) = ',nm,k); end
     info_1(A(k),s,cflag);
  end

end

% -------------------------------------------------------------------- %

function info_1(A,s,cflag)

  if isempty(A.Q) && isempty(A.data)
     fprintf(1,'\n%s (empty QSpace)',s);
     return
  end

  if ~isempty(s), fprintf(1,'\n%s\n',s);
  elseif ~cflag, fprintf(1,'\n'); end

  if isempty(A.Q)
     fprintf(1,'     Q:  []');
  else
     fprintf(1,'     Q:  %sx [%s]', ...
        vec2str(cellarrsize(A.Q,1),'-f'), ...
        vec2str(cellarrsize(A.Q,2),'-f') ...
     );
  end

  s={};
  if isfield(A.info,'otype') && ~isempty(A.info.otype)
     s{end+1}=A.info.otype;
  end
  if isfield(A.info,'itags') && ~isempty(A.info.itags)
     s{end+1}=[ itags2str(A.info.itags) ];
  end

  if ~isreal(A), s{end+1}='COMPLEX'; end

  if ~isempty(s), s(2,:)={',   '}; s=s([2 1],:); end
  s=cat(2,s{:});

  if ~isfield(A.info,'qtype')
       fprintf(1,'  (all abelian)%s\n',s);
  else fprintf(1,'  having ''%s''%s\n',A.info.qtype,s); end

  m=length(A.Q);
  s=A.data; s=whos('s'); s=num2str2(s.bytes,'-bytes');

  s={s,''};
% if isfield(A.info,'itags') && ~isempty(A.info.itags)
%      s{2}=['   ' itags2str(A.info.itags) ];
% else s{2}=''; end

  if isempty(A.data) || isempty(A.data{1})
     fprintf(1,'  data:  %d-D %s (%s) %s\n\n', m, class(A.data), s{:});
  else
     if ~isempty(A.Q) && ~isempty(A.Q{1})
        dd=getDimQS(A);
        if isvector(dd)
           dstr=vec2str(dd,'sep',' x ','fmt','%d','-f');
        else
         % dstr=mat2str2(dd,'sep',' x ','fmt','%d','rowsep',' => ','nofac');
           dstr=vec2str(dd(1,:),'sep',' x ','fmt','%d','nofac','-f');
           for i=2:size(dd,1)
               q=int2str2(dd(i,:)); dstr=[dstr, ' => ', strhcat(q,' x ')]; 
           end
        end
     else dstr=''; end

     fprintf(1,'  data:  %d-D %s (%s)      %s%s\n', ...
       m, class(A.data{1}), s{1}, dstr, s{2});
     if ~cflag, fprintf(1,'\n'); end
  end

end

% -------------------------------------------------------------------- %

