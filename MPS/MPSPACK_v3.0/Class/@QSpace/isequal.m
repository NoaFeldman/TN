function i=isequal(A,B)
% overloading isequal
% Wb,Nov09,09

  if isobject(A), A=struct(A); end
  if isobject(B), B=struct(B); end

  if ~isstruct(A) || ~isstruct(B) || numel(A)~=numel(B)
     i=0; return
  end

  if 0
   % for backward compatibility
   % NB! info field if empty can be either [] or ()
     if ~isempty(A) && ~isfield(A,'info'), A(1).info=[]; end
     if ~isempty(B) && ~isfield(B,'info'), B(1).info=[]; end

     for i=1:numel(A), if isequal(A(i).info,{}), A(1).info=[]; end, end
     for i=1:numel(B), if isequal(B(i).info,{}), B(1).info=[]; end, end
  end

  i=builtin('isequal',A,B);

end

