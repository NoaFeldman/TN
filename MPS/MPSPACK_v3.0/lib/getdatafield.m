function A = getdatafield(varargin)
% Function: getdatafield - get numerical data from field of structure array
% Usage: A = getdatafield(S, fieldname [,index])
%
%   if field elements actually contain arrays the optional
%   third argument index allows to specify which element to take.
%   index must be numeric except for the keyword 'end'
%
% Options
%
%   '-q'  quiet flag.
%
% Wb,May19,06

  getopt('init',varargin);
     qflag=getopt('-q');  % quiet (e.g. no NaN warning)

   % default value (if empty)
     if getopt('-0'), dval=0;
     else dval =getopt('dval',[]); end 

  varargin=getopt('get_remaining');

  if ~isempty(dval), qflag=1; end

  narg=length(varargin);
  if narg<2 || narg>3
     eval(['help ' mfilename]); A=[];
     if nargin, error('Wb:ERR','invalid usage'); else return; end
  end

  S=varargin{1};
  field=varargin{2};
  if narg<3, idx=[]; else idx=varargin{3}; end % index

  sa=size(S); e=0;
  na=prod(sa); if isempty(S), return; end

  if ~isfield(S,field), error('Wb:ERR', ...
     'invalid scalar field name (%s)',field); end

  f=getfield(S(1),field); isn=isnumeric(f);
  if isempty(idx) % && ~isscalar(f)
   % may have encountered empty entry
     for i=1:numel(S), if numel(getfield(S(i),field))>1
      error('Wb:ERR','index expected for data array'); end
     end
  end

  if isn 
     if isempty(dval), A=nan(sa); else A=repmat(dval,sa); end
  else A=repmat(f,sa); end

  C=cell(sa); eval(sprintf('C(:)={S.%s};',field));

  if isempty(idx)
     for i=1:na, if isempty(C{i}), e=1; break, end, end
     if e==0, A=reshape(catc(1,C{:}),sa);
     else, e=0;
      % structure array contains empty fields []
        for i=1:na
         % got empty or isscalar
           if numel(C{i})==1, A(i)=C{i}; else e=e+1; end
        end
     end
  else
     if ~ischar(idx)
        for i=1:na, try, A(i)=C{i}(idx); catch, e=e+1; end; end
     elseif isequal(idx,'end')
        for i=1:na, try, A(i)=C{i}(end); catch, e=e+1; end; end
     else
        wblog('ERR','Invalid index as arg #3'); idx
        return
     end
  end

  if e && ~qflag
     if isn, wblog('WRN','%d values set to NaN.',e);
     else wblog('WRN','%d structures are empty.',e); end
  end

end

