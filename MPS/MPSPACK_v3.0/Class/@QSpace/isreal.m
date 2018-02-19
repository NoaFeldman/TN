function isr=isreal(A)
% function i=isreal(A)
%
%    checks whether A contains complex data
%    (this overloads matlabs isreal() function).
%
% Wb,Apr23,15

% NB! use isreal rather than iscomplex, since isreal
% also corresponds to a standard matlab function.

  n=numel(A); isr=ones(size(A));

  if ~builtin('isempty',A)
     for k=1:n, dd=A(k).data;
        for i=1:numel(dd)
           if ~isreal(dd{i}), isr(k)=0; break; end
        end
     end
  end

end

