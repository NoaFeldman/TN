function display(A,varargin)
% Function display(A [,OPTS])
% Options
%
%   '-a'       show all entries
%   '-v'       verbose flag (shows CGS dimension separately for each symmetry)
%   '-E'       sort QSpace wrt. data (energy) [assumes diagonal Hamiltonian]
%   '-s'       sort QIDX
%   'sperm',.. sort QIDX using given permutation
%
% Wb,Mar01,08

  getopt('init',varargin);
     m  = getopt('m',8);
     nm = getopt('nm',''); % Wb,Apr08,14
     vflag=getopt('-v');
     aflag=getopt('-a'); if aflag, m=Inf; end % show all entries

     if     getopt('-E'), Eflag=1; os={'-E'};
     elseif getopt('-R'), Eflag=2; os={'-R'};
     elseif getopt('-S'), Eflag=3; os={'-S'};
     else Eflag=0;
         sflag=getopt('-s');
         sperm=getopt('sperm',{});
     end

  n=getopt('get_last','');

  if isempty(nm), nm=n;
  elseif ~isempty(n), {nm,n}
     error('Wb:ERR','invalid usage (name specified twice !?)')
  end

  if Eflag, [A,isd]=sort(A,os{:}); Eflag=isd;
  else
     if ~isempty(sperm)
        if  length(A(1).Q)~=length(sperm), sperm
            error('Wb:ERR','\n   ERR invalid sort-permutation'); end
        sflag=1; sperm={sperm};
     end
     if sflag, A=sort(A,sperm{:}); end
  end

  sA=size(A); n=prod(sA);
  if length(sA)==2 && sA(2)==1, sA=sA(1); end

  if isempty(nm)
     nm=inputname(1);
   % if isempty(nm) && vflag, nm='ans'; end
  end
  eflag=0; % empty flag

  if n==1
     if ~isempty(nm)
        if any(nm~=' ')
             s={sprintf('%s = ',nm)};
        else s={' '}; end
     else s={}; end
     display_1(A,m,Eflag,vflag,s{:});
  elseif n>1
     for i=1:n
        if i>1 && i<n && ise(A(i-1)) && ise(A(i+1))
           if ~eflag, fprintf(1,'\n...'); end
           eflag=eflag+1; continue;
        else eflag=0; end

        si=vec2str(ind2sub_aux(sA,i),'sep',',','-f');
        if ~isempty(nm)
             s=sprintf('%s(%s) = ',nm,si);
        else s=sprintf('%s>',si); end

        if ~ise(A(i))
           if n<8 || vflag % || i==1 || i==n
                display_1(A(i),m,Eflag,vflag,s);
           else info(A(i),s, 'c'); end
        else fprintf(1,'\n%s(empty QSpace)',s); end
     end
  else
     fprintf(1,'\n%s = \n', nm);
     fprintf(1,'  (empty QSpace)\n');
  end
  
  if n>=1 && ise(A(end))
  fprintf(1,'\n\n'); else fprintf(1,'\n'); end

end

function i=ise(A)
  % i = isempty(A) || isempty(A.Q) || isempty(A.Q{1});
    i = isempty(A) || isempty(A.data); % be aware of scalars!
end

% -------------------------------------------------------------------- %
% -------------------------------------------------------------------- %
% calculate index zero based first, then convert
% NB! ind2sub() would not return vector, but variable number
% of output arguments ranging from 1 to length of s

function kk = ind2sub_aux(s,k)

 % if vector keep it as vector (i.e. single index)
   if numel(s)<=2 && sum(s>1)<2, kk=k; return; end

   k=k-1; kk=zeros(size(s));

   for i=1:length(s)
   kk(i)=mod(k,s(i)); k=(k-kk(i))/s(i); end

   kk=kk+1;

end

% -------------------------------------------------------------------- %
% -------------------------------------------------------------------- %
% introduce compact notation for SU(N), Sp(2n), etc
% Wb,Apr02,15

% see also QSpace::getqfmt() // Wb,Dec14,15

function qfmt=get_q_fmt(qtype,r,m)

   q=strread(qtype,'%s','whitespace',' ,;|\n\r\t')';
   if isempty(q)
      qfmt=repmat(' %2g',1,m);
      if r>1, qfmt=repmat([ ' ;' qfmt],1,r); end
      qfmt=qfmt(4:end);
      return
   end

   for i=1:numel(q)
      j=regexp(q{1,i},'\d');
      if isempty(j), q{2,i}='%2g';
      else
         w=q{1,i}(1:j-1); n=str2num(q{1,i}(j:end));
         switch w
            case 'SU', q{2,i}=repmat('%g',1,n-1);
            case 'Sp', q{2,i}=repmat('%g',1,n/2);
            case {'Z','P'}, q{2,i}='%2g';
            otherwise
            error('Wb:ERR','\n   ERR unexpected symmetry %s',q{1,i});
         end
      end
   end
   q(3,:)={' '}; q{3,end}='';

   qfmt=[q{2:3,:}];

   if nargin>1
      qfmt=repmat({qfmt},2,r);
      qfmt(2,:)={' ; '}; qfmt{end}='';
      qfmt=[qfmt{:}];
   end

end

% -------------------------------------------------------------------- %
% -------------------------------------------------------------------- %

function display_1(A,m,Eflag,vflag,varargin)

  info(A,varargin{:}); r=length(A.Q);
  dlen=numel(A.data);

  cgflag=gotCGS(A);

% size format
  if cgflag
     ns=length(find(A.info.qtype==','))+1;
     if ~isempty(A.Q) && ~isequal(size(A.info.cgr),[dlen, ns])
     error('Wb:ERR','CG size mismatch'); end
     nq=size(A.info.cgr,2);
     sfmt=sprintf('%%-%ds', 3*length(A.Q));
  else
     sfmt=sprintf('%%-%ds', 5+4*length(A.Q));
  end

  NL=1E10; % deactivate

  fstr={ '%7.4g'       % for scalars (allow for complex numbers)
         '  %-12s%s%s\n' };  % for non-scalars: size in bytes

  if isequal(get(0,'Format'),'long')
  fstr{1}='%16.13g'; end

  if ~isempty(A.Q)
     QQ=permute(cat(3,A.Q{:}),[3 2 1]);
     if ~isempty(A.info), t=A.info.qtype; else t=''; end
     qfmt=get_q_fmt(t,numel(A.Q),size(A.Q{1},2));
  else QQ=[]; end

  i=0; rflag=1; % m=8;
  while i<dlen, i=i+1;
      if rflag && dlen>20 && i>=m
         fprintf(1,'\n  (skipping %d entries)\n\n',dlen-2*m+1);
         i=dlen-m; rflag=0;
         continue
      end

      if i>dlen-min(2,NL) || mod(i,NL), nl=''; else nl='\n'; end
      Ai=A.data{i}; sa=size(Ai); % sa(end+1:m)=1;
      if cgflag
         sc=cell(1,nq);
         for j=1:nq, sc{j}=cgr_size(A,i,j); end
         sc=cat2(1,sc{:},{1}); sc(:,end+1:r)=1; sa(end+1:r)=1;

         s1=vec2str(sa,'fmt','%g','sep','x','-f');
         if ~vflag
            sc=prod(sc,1);
            s2=vec2str(sc,'fmt','%g','sep','x','-f');
            s2=sprintf(sfmt,s2);
         else
            n=size(sc,1); s2=cell(1,n);
            for j=1:n
               s2{j}=vec2str(sc(j,:),'fmt','%g','sep','x','-f');
            end
            s2=sprintf(' %6s',s2{:}); s2=s2(2:end);
         end

         fprintf(1,['%6d.  ' sfmt ' | %s' ],i,s1,s2);
      else
         s=vec2str(sa,'fmt','%g','sep','x','-f');
         fprintf(1,sfmt, sprintf('%6d.  %s',i,s));
      end

      if ~isempty(QQ)
         % s=mat2str2(QQ(:,:,i),'rowsep',' ; ','fmt','%2g','-f');
           s=sprintf(qfmt,QQ(:,:,i)');
      else s=''; end
      fprintf(1,' [ %s ]',s);

      s=prod(sa); dfac=1;
      if 1 || Eflag %  || s==1
         if cgflag
          % add CG range-factor only if != 1
            [dfac,sc]=get_cgr_fac(A,i); % cgw
            if ~isempty(sc); sc=['{' sc '}']; end
         else sc=''; end
      end

      if s==1, Ai=dfac*Ai;
         if isreal(Ai)
              q=sprintf(fstr{1}, Ai);
              if isempty(find(q=='.' | q=='e',1)), q=[q '.']; end
              fprintf(1,[';  ' q '    %s\n' nl], sc);
         else fprintf(1,[';  %s %s\n' nl], num2str2(Ai,'fmt',fstr{1}), sc);
         end
         continue

      elseif Eflag
         if ~isreal(Ai), error('Wb:ERR',...
           '\n   ERR real data{} expected (got complex)'); end
         if Eflag==1, Ai=diag(Ai); sa=size(Ai); end
         if numel(find(sa>1))>1, error('Wb:ERR',...
           '\n   ERR invalid data{} (diagonal representation expected)'); end
         n=numel(Ai);
         if n<4
            q=sprintf([' ' fstr{1}],Ai);
         else q=[ ...
            sprintf([' ' fstr{1}],Ai(1:2)), sprintf(' ..(%d).. ',n-3), ...
            sprintf(fstr{1},Ai(end)) ];
         end
         fprintf(1,[';  %s %s\n' nl], q(2:end), sc);
         continue
      end

      s=s*8; % assume double size

      if     s<2^10, s=sprintf('%g B',s);
      elseif s<2^20, s=sprintf('%.1f kB',s/2^10);
      else           s=sprintf('%.1f MB',s/2^20);
      end

      s={s,'',''};
      if ~isempty(sc), s{2}=[sc]; end
    % NB! add 'COMPLEX' info to header as printed by info() above!
    % Wb,Apr23,15
    % if ~isreal(A.data{i}), s{3}='  complex'; end

      fprintf(1,[fstr{2} nl],s{:});
  end

% fprintf(1,'\n')

end

% -------------------------------------------------------------------- %
% -------------------------------------------------------------------- %
