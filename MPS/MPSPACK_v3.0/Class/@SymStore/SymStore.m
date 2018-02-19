function q=SymStore(sym,task,varargin)
% SymStore is a matlab class to analyze the RCStore for a specfic symmetry
% Wb,Jan30,15

  sym=regexprep(sym,'(S[pU])\(*(\d+)\)','$1$2');
  D=get_dir(sym); % just double check whether RCStore exists

  if nargin<2 || ~ischar(task)
     helpthis, if nargin || nargout
     error('Wb:ERR','invalid usage'), end, return
  end

  if isequal(task,'dim')
     getSymDim(sym,varargin{:});
  elseif isequal(task,'XOM')
     getSymXOM(sym,varargin{:});
  elseif isequal(task,'loadC')
     q=LoadCData(sym,varargin{:});
  else error('Wb:ERR','\n   ERR invalid task ''%s''',task);
  end

