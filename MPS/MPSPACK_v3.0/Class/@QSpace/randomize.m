function A=randomize(A)
% function A=randomize(A)
% Wb,Nov11,09

  n0=0; n2=0; n=numel(A.data); 

  for i=1:n, d=A.data{i};
     n0=n0+sum(reshape(d.*d,[],1)); d=randn(size(d));
     n2=n2+sum(reshape(d.*d,[],1)); A.data{i}=d;
  end

  nfac=sqrt(n0/n2);
  for i=1:n, A.data{i}=A.data{i}*nfac; end

end

