function C=plus(A,B)
% overloading the + operator

  if isnumeric(B) && isscalar(B) && length(A.Q)==2
   % adding single scalar = acts like diagonal operator
     [isd,s]=isdiag(A);
     if ~isempty(s), error('Wb:ERR',estr); end

     C=A; n=length(C.data);

     if isd>1
      % assume compressed representation
        for i=1:n, C.data{i} = C.data{i} + B; end
     else
        for i=1:n
           if isequal(C.Q{1}(i,:),C.Q{2}(i,:))
              C.data{i} = C.data{i} + full(B*speye(size(C.data{i})));
           end
        end
     end
  else
     if ~isequal(size(A),size(B)), error('Wb:ERR',...
        '\n   ERR invalid usage (size mismatch between A and B)'); end
     C=A; n=numel(A);
   % QSpace() required by matlab/2016a // Wb,Aug03,16
     for i=1:n, C(i)=QSpace(plusQS(A(i),B(i))); end
  end

% some routines require all blocks including zeros (eg. setup routines)
% C(1)=skipZerosQS(C);

end

