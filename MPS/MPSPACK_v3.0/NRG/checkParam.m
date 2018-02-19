function e=checkParam(strEval)
% checkParam(strEval) - update parameters to struct()
%
% WRN
%
%    strEval evaluates string expressions (default: 1)
%    this might result in circular function calls
%    -> need to unset.
%
% Wb,Oct11,05

  global param paras

  e=0; if nargin==0, strEval=1; end

  if ~isempty(param)
     if isstruct(param)
        if strEval
           ff=fieldnames(param);
           for i=1:length(ff)
               eval(['fval=param.' ff{i} ';']);

               % isempty(findstr(fval,'weichsel')) && ...
               if ~ischar(fval) || ~isempty(findstr(fval,'//'))
               continue; end

             % evaluate string (assuming it to be arithmetic)
               try eval([ 'param.' ff{i} '=' fval ';']);
               catch
                  if ~isempty(regexp(fval,'[*+-()]'))
                      e=e+1; wblog('ERR',...
                     'Failed to evaluate field string `%s''',fval);
                  end
               end
           end
        end
     elseif iscell  (param), param=param2struct(param);
     elseif isvector(param) & isreal(param)
        if length(param)>4, B=param(5); else B=0; end
        param=struct(...
             'Gamma',  param(1), ...
             'U',      param(2), ...
             'epsd',   param(3), ...
             'Lambda', param(4), ...
             'Bfield', B         ...
        );
     else
        istr='Invalid/unknown param structure!';
        if nargout, e=1; wblog('ERR',istr); return;
        else error('Wb:ERR',str); end
     end
  elseif ~isempty(paras)
     param=struct(...
          'Gamma',  paras(3), ...
          'U',      paras(2), ...
          'epsd',   paras(1), ...
          'Lambda', paras(4)  ...
     );
  else
     istr='Neither paras nor param are set global';
     if nargout, e=1; wblog('ERR',istr); return;
     else error('Wb:ERR',istr); end
  end

  if ~nargout, clear e; end

end

