function C=mtimes(A,B,ia,ib)
% function C=mtimes(A,B,ia,ib)
% overloading * operator

  if isnumeric(A) || isnumeric(B)

     if nargin~=2, error('Wb:ERR','Invalid number of arguments'); end
     if     isnumeric(A), fac=A; C=B;
     elseif isnumeric(B), fac=B; C=A; end

     for k=1:numel(C), n=length(C(k).data);
         data=C(k).data; for i=1:n, data{i}=fac*data{i}; end
         C(k).data=data;
     end

  else
   % similar usage than contractmat(A,M,ia)

     if nargin<2, error('Wb:ERR','Invalid number of arguments'); end
     if isempty(A.Q) && isempty(A.data) C=QSpace(); return, end
     if isempty(B.Q) && isempty(B.data) C=QSpace(); return, end

     ra=length(A.Q); if nargin<3, ia=[]; end % plain matrix product A*B
     rb=length(B.Q); if nargin<4, ib=[]; end % plain matrix product A*B

     if nargin<3 && ra==3 && rb==3 % Wb,Dec04,15
        qa=uniquerows(A.Q{3}); qb=uniquerows(B.Q{3});
        if size(qa,1)==1 && size(qb,1)==1
           if ~isequal(qa,qb), error('Wb:ERR',...
              '\n   ERR invalid usage (operator mismatch !?)'); 
           end
           ia='13*'; ib='13';
        end
     end

     if isempty(ia) % plain matrix product A*B
        if ra<2 || ra>3 || ~isop(A), error('Wb:ERR',...
           '\n   ERR invalid usage (A must be operator)'); end
        ia=2;
     end
     if isempty(ib) % plain matrix product A*B
        if rb<2 || rb>3 || ~isop(B), error('Wb:ERR',...
           '\n   ERR invalid usage (B must be operator)'); end
        ib=1;

        if ra>2 && rb>2, error('Wb:ERR',...
           '\n   ERR invalid usage (got two rank-3 IROPs) !?');
         % presumably need to contract A'*B (i.e. including dagger!)
         % Wb,Aug27,15
        end
     end

     if ~ischar(ia) && ia>ra || ~ischar(ib) && ib>rb, error('Wb:ERR',...
        'index out of bounds (A: %d/%d; B: %d/%d)',ia,ra,ib,rb); end

     if ra==3 && rb==2
        % p=[ 1:ia-1, ra, ia:ra-1 ]; // Wb,Aug27,15
          p=[ 1 3 2 ];
     else p=[]; end

     C=contractQS(A,ia,B,ib,p);
     if numel(C.Q)==3, C.info.otype='operator'; end

     C=class(C,'QSpace');

  end

end

