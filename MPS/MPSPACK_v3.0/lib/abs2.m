function d2 = abs2(dd)
% Function: d2 = abs2(dd)
%
%    calculate |dd|� (elementwise)
%
% Wb,Apr30,04

  d2 = dd .* conj(dd);

end

