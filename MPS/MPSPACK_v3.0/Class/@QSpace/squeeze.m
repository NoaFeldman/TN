function [A,I]=squeeze(A,varargin)
% Function [A,I]=squeeze(A [,OPTS])
%
%    remove singleton dimensions (similar to matlab's routine)
%
% Options
%
%   '-start'  only singletons at the front
%   '-end'    only singletons at the end (default: all)
%
% Wb,Apr24,08

  if nargin>1 && isnumeric(varargin{1})
     ks=varargin{1}; varargin=varargin(2:end);
     narg=length(varargin);
  else ks=[]; end

  getopt('init',varargin);
     sflag=getopt('-start');
     eflag=getopt('-end'  );
  getopt('check_error');

  if isempty(ks)
   % by default skip singleton(s) at end of QSpace
     if ~sflag && ~eflag, eflag=1; end
  elseif sflag || eflag, wblog('WRN',...
    '-start and -end ignored since index is explicitely specified');
  end

  m=zeros(size(A));
  for k=1:numel(A)
     Ak=struct(A(k)); s=dim(A(k),'-f'); r=numel(Ak.Q); is=[]; p=[];

     if ~isempty(ks)
        if any(ks>r) || any(s(ks)~=1), error('Wb:ERR',['given index ' ...
           '[%s] does not correspond to singleton'],vec2str(ks)); end
        is=ks;
     else
        if sflag
           for i=1:r, if s(i)~=1, is=1:(i-1); break; end, end
        elseif eflag
           for i=r:-1:1, if s(i)~=1, is=(i+1):r; break; end, end
        end
     end

     if isempty(is)
        if nargout>1
           I(k)=add2struct('-',s,is,p);
        end
        continue
     end

     Ak.Q(is)=[];
     
   % permute singletons to end
     p=1:r; pk=p(is); p(is)=[]; p=[p,pk];

     n=numel(Ak.data); for i=1:n
     Ak.data{i}=permute(Ak.data{i},p); end

     if isfield(Ak.info,'cgs') && ~isempty(Ak.info.cgs)
        if isstruct(Ak.info.cgs), cgs=Ak.info.cgs;
           l=max(2,r-numel(is));
           ik=p(1:l); i2=p(l+1:end); n=numel(cgs);
           for i=1:n
              if any(cgs(i).S(i2)~=1), error('Wb:ERR',...
                '\n   ERR got CGS dimension!=1 !??'); end
              cgs(i).S=cgs(i).S(ik);
           end
           Ak.info.cgs=cgs;
      % else sparse array already has minimum rank of 2!
        end
     end
     A(k)=Ak;

     if nargout>1
        I(k)=add2struct('-',s,is,p);
     end
  end

end

