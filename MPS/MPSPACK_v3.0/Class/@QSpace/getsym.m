function [ss,s]=getsym(A,varargin)
% function [ss,s]=getsym(A [,dim,opts])
% Options
%
%   -c        return symmetries split as cell array
%   -d        get dimensions (rank) of each symmetry
%   'ir',dim  index range in q-labels for given dimension
%   '-I',dim  get extended info structure for given symmetry
%
% Wb,Oct24,11

  iflag=0;
  getopt('init',varargin);
     cflag=getopt('-c');
     dflag=getopt('-d');

     ir=getopt('ir',[]); if isempty(ir)
     ir=getopt('-I',[]); if ~isempty(ir), iflag=1; end; end % Wb,Apr11,12

  dim=getopt('get_last',[]);

  if ~isempty(dim)
     if ~isnumber(dim) || cflag
        wblog('ERR','%s() invalid usage',mfilename);
        return
     end
     cflag=1;
  end
  if dflag || ~isempty(ir), cflag=1; dflag=1; end

  if dflag, ss=[]; elseif cflag, ss={}; else ss=''; end
  if isempty(A), return; end

  for k=1:numel(A), Ak=struct(A(k));
     if isfield(Ak.info,'qtype') && ~isempty(Ak.info.qtype), s=Ak.info.qtype;
        q=strread(s,'%s','whitespace',', '); % safeguard
        if numel(q)>size(Ak.Q{1},2), error('Wb:ERR',['\n   ERR ' ...
       'inconsistent symmetry type (''%s''; got rank-%d)'],s,size(Ak.Q{1},2)); end
        if cflag, s=q; end
     else
        if ~isempty(Ak.Q), ns=size(Ak.Q{1},2);
           if cflag, s=repmat({'A'},1,ns);
           else s=repmat(',A',1,ns); s=s(2:end); end
        else s={}; end
     end
     if k>1, if ~isequal(s,ss), error('Wb:ERR',...
       '\n   ERR inconsistent symmetry type (%s <> %s)',s,ss); end
     else ss=s; end
  end

  if dflag
     for i=1:numel(ss), s=ss{i};
        if isequal(s,'A'), ds(i)=1;
        elseif ~isempty(regexp(s,'^SU\d+$')), ds(i)=str2num(s(3:end))-1;
        elseif ~isempty(regexp(s,'^Sp\d+$')), ds(i)=str2num(s(3:end))/2;
        else error('Wb:ERR','\n   ERR unknown symmetry `%s''',s);
        end

        if ds(i)~=round(ds(i)), error('Wb:ERR',...
           '\n   ERR invalid symmetry rank %g (%s)',ds(i),s);
        end
     end

     if ~isempty(ir)
        if numel(ir)>1, error('Wb:ERR',['\n   ' ... 
          'ERR SINGLE index expected for ''ir''']); end
        if ~isempty(dim), error('Wb:ERR',[ '\n   ' ...
          'ERR dimension shall be specified with options `ir'' itself']); end

        s=ss{ir};
        ic=cumsum([1 ds]); ss=ic(ir):(ic(ir+1)-1);

        if iflag
         % is: which symmetry (i.e. 1st, 2nd, ...)
         % i : corresponding column indizes within Q(:,i)
         % r : rank of symmetry
           ss=struct('sym',s,'is',ir,'j',ss,'r',numel(ss));
           clear s
        end
     else
        if isempty(dim), s=ss; ss=ds; else s=ss{dim}; ss=ds(dim); end
     end
  else
     if ~isempty(dim), ss=ss{dim}; end
     clear s
  end

end

