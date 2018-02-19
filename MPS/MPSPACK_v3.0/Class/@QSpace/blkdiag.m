function M=blkdiag(A,varargin)
% function M=blkdiag(A,varargin)
% Wb,Aug10,16

  if numel(A)~=1, error('Wb:ERR','invalid usage'), end

  if ~isempty(A.Q) % NB! A may also be a scalar
     if numel(A.Q)~=2 || ~isequal(sort(A.Q{1},2),sort(A.Q{2},2))
      % quick and dirty check for scalar operators
      % or wave functions with all-in labels (strictly would
      % have to check for dual irreps)
        wblog('WRN','got rank-%g object with off-diagonal blocks !?',...
        numel(A.Q));
     end
  end

  M=A.data; M=blkdiag(M{:});
  
end

