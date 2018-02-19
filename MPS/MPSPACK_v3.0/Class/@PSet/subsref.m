function [p,pstr,tdone,jj]=subsref(P,SS)
% Function [p,pstr,tdone,jj]=subsref(P,S)
% Overloading the subsref() routine for PSet. Usage:
%
%  * Regular call by serial index mean for 
%       for ip=1:P.n, [p,pstr]=P(ip); ...; end
%
%  * Add time finished estimate in parameter loop: add flag -time
%       for ip=1:P.n, [p,pstr]=P(ip,'-time'); ...; end
%       for ip=1:P.n, [p,pstr,tdone]=P(ip); ...; end
%    
%  * At an intermediate stage reference to a specific parameter
%    can be made by its name; this will give the index iz and its
%    maximum value nz for parameter 'B'
%
%       [iz,nz]=P(ip,'B');
%
% Output
%    
%    p     specific parameters setting for given index ip
%          returned as structure (array, if length(ip)>1)
%    pstr  info string that can be used for logging (for scalar ip only).
%    tdone estimated time that job will be finished (extrapolate ip->P.n)
%    jj    specific index set (like ind2sub)
%
% Wb,Jul31,08

  S=SS(1); ns=length(SS);

  if isequal(S.type,'()'), sflag=0; n=length(S.subs);
     if n, e=0;
        if n==2, if ischar(S.subs{2}), sflag=1; else e=1; end, end
        if e || n>2 || ~isnumeric(S.subs{1})
        error('Wb:ERR','invalid usage'); end
     else S.subs{1}=P.i; end

     ii=reshape(S.subs{1},1,[]); ni=length(ii); jj=zeros(P.r,ni);

     if any(ii>P.n) || isempty(ii), error('Wb:ERR',...
     'index out of range (%g/%g)',max(ii),P.n); end

     if sflag && isequal(S.subs{2},'-time'), sflag=-1; end


     i0=ii;

     jj=cell(1,numel(P.s)); [jj{:}]=ind2sub(P.s,i0);
     jj=cat(2,jj{:});

     if sflag>0, vn=S.subs{2};
        for i=1:P.r, if isequal(P.name{i},vn)
            p=jj(:,i); pstr=P.s(i); return;
        end,end
        p=0; pstr=0; return
     end

     s=P.name;
     for i=size(jj,1):-1:1
         for j=1:P.r, s{2,j}=P.data{j}(jj(i,j)); end
         p(i)=struct(s{:});
     end

     if ns>1, p=builtin('subsref',p,SS(2:end)); end

     if nargout>1
        if ni>1, wblog('WRN',...
        'pstr only specified for SINGLE input index (%g)',ni); end

        s(4,:)=s(2,:);
        s(2,:)=mat2cell(jj(1,:),1,ones(1,P.r));
        s(3,:)=mat2cell(P.s,    1,ones(1,P.r));

        for i=1:P.r
           if s{3,i}>1
                s{1,i}=sprintf('%s(%g/%g)=%.4g  ', s{:,i});
           else s{1,i}=sprintf('%s=%.4g  ', s{[1 4],i}); end
        end

        pstr=[s{1,:}]; pstr=pstr(1:end-2);
        
        if numel(i0)~=1, return; end

        if sum(P.s>1)>1, pstr=[ sprintf('%2g/%g:  ',i0,P.n), pstr ]; end

        if sflag>=0 && nargout<3, return; end

        P.t(end+1,:)=[i0, now]; n=size(P.t,1);
        if P.s(1)<=1, m=11; else; m=P.s(1)*ceil(20/P.s(1))-1;
        end
        P.t=P.t(max(1,n-m):n,:);

        if n>1 && all(diff(P.t(:,1))==1)
           if size(P.t,1)>6
             np=2;

             [pp,s,mu]=polyfit(P.t(2:end,1),log(diff(P.t(:,2))),np);
             [t,dt]=polyval(pp,(i0:P.n),s,mu);
             t=exp(t);
             if all(abs(dt)<1)
                  dt=24*60*(t*exp(dt)'); t=now+sum(t);
             else dt=[]; end
           else dt=[]; end
           if isempty(dt)
              m=(P.n-i0+1);
              t=P.t(end,2)+m*mean(diff(P.t(:,2)));
              dt=sqrt(m)*24*60*std(diff(P.t(:,2)));
           end

           if dt>1
              if dt>99, dt=dt/60;
                   s=sprintf(' (@ %.2g hrs)',dt);
              else s=sprintf(' (@ %.2g min)',dt); end
           else s=''; end

           tdone=sprintf('estimated time finished: %s%s',datestr(t),s);
           if nargout<3, pstr=[pstr, char(10), tdone];  end

        else tdone=''; end

        n=inputname(1); if ~isempty(n), assignin('caller',n,P); end
     end
     
  elseif isequal(S.type,'.') && ischar(S.subs)

     if builtin('isfield',struct(P),S.subs)
        p=builtin('subsref',P, S);
        if ns>1, p=builtin('subsref',p,SS(2:end)); end
        return
     end

     n=S.subs;
     for i=1:P.r
        if isequal(P.name{i},n), p=P.data{i};
           if ns>1, p=builtin('subsref',p,SS(2:end)); end
           return
        end
     end

     error('Wb:ERR','name ''%s'' not declared within PSet',S.subs); 

  else
     p=builtin('subsref',P,SS);
  end

end

