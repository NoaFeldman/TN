function C=comm(A,varargin)
% function C=comm(A [,B,cflag])
%
%    cflag if set, must equal 'conj'
%    which is interpreted as conjA, i.e. C = {A',B}
%
%    For anticommutator, use comm(A [, B, cflag], '+');
%
%    Here the operators must be of rank-2 or rank-3;
%    in the case of rank-3, this assumes the index order (ss',o),
%    i.e. the irop (=spinor) index comes last.
%
% Wb,Feb06,16

  if nargin>1 && ~ischar(varargin{1})
     B=varargin{1}; [conjA,wc]=get_opts(varargin{2:end});
     na=numel(A); nb=numel(B); C=QSpace(na,nb);
     for i=1:na
        for j=1:nb, C(i,j)=commAB_1(A(i),B(j),conjA,wc); end
     end

     if     na==1, C=reshape(C,size(B));
     elseif nb==1, C=reshape(C,size(A));
     end
  else
     C=QSpace(size(A)); na=numel(A);
     [conjA,wc]=get_opts(varargin{:});
     for i=1:na, C(i)=commAA_1(A(i),wc); end
  end

end

% -------------------------------------------------------------------- %

function [conjA,wc]=get_opts(varargin)

  conjA=0; wc=-1; if ~nargin, return; end

  getopt('init',varargin);
     conjA=getopt('conj');
  q=getopt('get_last',[]);

  if ~isempty(q) 
     if ~isequal(q,'+'), error('Wb:ERR',...
       '\n   ERR invalid usage (expecting ''+'' for anticommutator)'); end
     wc=+1;
  end

end

% -------------------------------------------------------------------- %
 
function C=commAB_1(A,B,conjA,wc)

  ra=numel(A.Q);
  rb=numel(B.Q);

  if ra<2 || ra>3, error('Wb:ERR',...
     '\n   ERR invalid usage (got A of rank %g !?)',ra); end
  if rb<2 || rb>3, error('Wb:ERR',...
     '\n   ERR invalid usage (got B of rank %g !?)',rb); end

  if ra==rb
     if conjA
        if ra==2
             C=[contractQS(A,'1*', B,'1' ), contractQS(B,'2', A,'2*' )];
        else C=[contractQS(A,'13*',B,'13'), contractQS(B,'23',A,'23*')];
        end
     else
        if ra==2
             C=[contractQS(A,'2', B,'1' ), contractQS(B,'2', A,'1' )];
        else C=[contractQS(A,'23',B,'13'), contractQS(B,'23',A,'13')];
        end
     end
  else
     if conjA
          C=[contractQS(A,'1*',B,'1'), contractQS(B,'2',A,'2*')];
     else C=[contractQS(A,'2', B,'1'), contractQS(B,'2',A,'1' )];
     end
     if ra==2, i=2; else i=1; end
     C(i)=permuteQS(C(i),'132');
  end

  if wc<0
       C=QSpace(C(1))-C(2); % regular commutator
  else C=QSpace(C(1))+C(2); % (fermionic) anticommutator
  end

  if ra~=rb
   % got rc=3 => inherit otype from A or B (if set)
     if ra>2, q=A.info.otype; else q=B.info.otype; end
     C.info.otype=q;
  end

end

% -------------------------------------------------------------------- %

function C=commAA_1(A,wc)

  ra=numel(A.Q);
  if ra<2 || ra>3, error('Wb:ERR',...
    '\n   ERR invalid usage (got A of rank %g !?)',ra); end
  
% always uses conjA
  if ra==2
       C=[contractQS(A,'1*', A,'1' ), contractQS(A,'2', A,'2*' )];
  else C=[contractQS(A,'13*',A,'13'), contractQS(A,'23',A,'23*')];
  end

  if wc<0
       C=QSpace(C(1))-C(2); % regular commutator
  else C=QSpace(C(1))+C(2); % (fermionic) anticommutator
  end

end

% -------------------------------------------------------------------- %

