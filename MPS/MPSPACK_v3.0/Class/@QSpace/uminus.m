function A=uminus(A)
% overloading the + operator

  for k=1:prod(size(A))
     for i=1:length(A(k).data)
        A(k).data{i}=-A(k).data{i};
     end
  end

return

