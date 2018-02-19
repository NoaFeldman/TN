function A=cellarrsum(C)
% Function: cellarrsum
% Usage: A=cellarrsum(C);
%
%    A = sum_i C{i}
%
% Wb,Nov30,05

  A=C{1}; n=numel(C);
  for i=2:n, A=A+C{i}; end

% turns out slower
% n=length(C(:)); s=size(C{1});
% A=reshape(cell2mat(C(:)'),[s n]);
% A=permute(sum(permute(A,[3 1 2])), [2 3 1]);

end

