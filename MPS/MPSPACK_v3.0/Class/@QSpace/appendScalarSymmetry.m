function A=appendScalarSymmetry(A,qtype,varargin)
% function A=appendScalarSymmetry(A,qtype [,opts])
% Wb,Jul12,11 ; Wb,Feb06,14 ; Wb,Dec08,14

% removed q from original usage [A=appendScalarSymmetry(A,q,qtype)]
% and added pos as optional argument // Wb,Feb06,14


  getopt('init',varargin);
     q  =getopt('q',[]);
     pos=getopt('pos',0);
  getopt('check_error');

  if isempty(q)
     if isequal(qtype,'A'), r=1;
     elseif ~isempty(regexp(qtype,'SU\d$'))
        r=str2num(qtype(end))-1;
     elseif ~isempty(regexp(qtype,'Sp\d$'))
        r=str2num(qtype(end))/2;
     else
       error('Wb:ERR','\n   ERR invalid/unknown symmetry (%s)',qtype);
     end
     if r<1 || r>9 || r~=round(r)
        error('Wb:ERR','\n   ERR invalid usage (got rank=%g !??)',r);
     end
     q=zeros(1,r);
  end

  for k=1:numel(A), Ak=A(k); Ik=Ak.info;
     qd_=getqdir(Ak,'-s');

     if ~isempty(Ik.cgr)
        if ~isstruct(Ik.cgr)
           error('Wb:ERR','\n   ERR invalid info.cgr field'); end
        Ik.qtype=[Ik.qtype ',' qtype];

        qd=Ik.cgr(1).qdir;
        if isempty(qd), qd=qd_; end

      % by keeping remaining fields (cid, size, cgt, ...) empty,
      % this either copies this information from the existing
      % corresponding CData in gCG.BUF if present, otherwise the
      % corresponding CData is loaded from the RCStore
      % see also QSpace/getvac // INIT_SCALAR // Wb,Jan29,15
        c=emptystruct(Ik.cgr);
          c.type=qtype;
          c.qset=repmat(q,1,numel(Ak.Q));
          c.qdir=qd;
          c.cgw=['1' 0];
        Ik.cgr(:,end+1)=c;
        Ak.info=Ik;

     elseif isequal(qtype,'A')
        if ~isempty(Ik.qtype)
           Ak.info.qtype=[Ak.info.qtype ',' qtype];
        end
     else
        if isempty(Ik.qtype)
             Ik.qtype=[repmat('A,',1,size(Ak.Q{1},2)), qtype];
        else Ik.qtype=[Ik.qtype ',' qtype];
        end
        qt=strread(Ik.qtype,'%s','delimiter',',;');
        nq=numel(qt);

        c=repmat(getRC('--empty'),1,nq);
          c(end).type=qtype;
          c(end).qset=repmat(q,1,numel(Ak.Q));
          c(end).qdir=qd_;
          c(end).cgw=['1' 0];
        Ik.cgr=repmat(c,size(Ak.Q{1},1),1);
        Ak.info=Ik;

     end

     Q=Ak.Q;
     for i=1:numel(Q), Q{i}=[Q{i},repmat(q,size(Q{i},1),1)]; end
     Ak.Q=Q;

     if pos>0, n=numsym(Ak);
        if pos>n, error('Wb:ERR',...
           '\n   ERR index pos out of range (%g/%g) !?',pos,n);
        elseif pos<n
           perm=[1:pos-1,n,pos:(n-1)];
           Ak=symperm(Ak,perm);
        end
     end

     A(k)=Ak;

  end

end

