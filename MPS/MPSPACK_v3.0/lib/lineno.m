function [line, name] = lineno(varargin)
% Function: lineno - get line-number and name of calling function
% Usage: [line, name] = lineno([id=1, OPTS])
%
%   returns line number and function name of calling function,
%   (0 and '' if called from the Matlab command-line prompt.)
%
%   id sets the offset - id=1 = current mfile, id=2 = calling mfile, ...
%
% Options
%
%   'nx',..  number of highest stack entries to skip
%   '-s'     short (do not print line numbers)
%   'all'    display full stack inline
%   'ALL'    display full stack inline (including subroutine names)
%   'len',.. max. length of output string (for nargout==1 only)
%   '<FMT>'  explicit format specification
%            '%.20F -> %S:%L'  (F=file, S=subroutine, L=line)
%
% Output - if nargout
%
%   =0: dbstack is displayed only
%   =1: returned object is single string "name:line"
%       an optional length argument shortens the string using shortfstr()
%
% Ref. http://woodshole.er.usgs.gov/...
% staffpages/cdenham/public_html/snackbar/lineno.m
%
% Wb,Nov13,05  Wb,Jun16,07

  if nargin && isnumeric(varargin{1})
  id =varargin{1}; varargin=varargin(2:end); else id=1; end

  getopt('INIT',varargin);
     if getopt('ALL'), aflag=2; 
     elseif getopt('all'); aflag=1; else aflag=0; end
     nx   =getopt('nx',0)+2; % skip *this routine
     sflag=getopt('-s');
     len  =getopt('len',-1);
  varargin=getopt('get_remaining');

  if length(varargin)+aflag>1
  error('Wb:ERR','invalid usage'); end

  if length(varargin), fmt=varargin{1}; else fmt=[]; end

  [stack, index] = dbstack;
  % stack = trace information
  % index = current workspace index

  if aflag
     if aflag==2
        line={stack.file; stack.name; stack.line};
        line=fliplr(line(:,nx:end));
        if sflag, line=line(1:2,:);
             line=sprintf('%s::%s -> ',line{:});
        else line=sprintf('%s::%s:%d -> ',line{:});
        end; line=line(1:end-4);
     elseif aflag==1
        line={stack.file; stack.line}; line=fliplr(line(:,nx:end));
        if sflag, line=line(1,:);
             line=sprintf('%s -> ',line{:});
        else line=sprintf('%s:%d -> ',line{:});
        end; line=line(1:end-4);
     end
     return
  end

  if length(stack) > 1
     id = max([min([1+id, length(stack)]), 1]);
     stack = stack(id);
  else
     stack=struct('line',0, 'file','', 'name','');
  end

  if ~isempty(fmt)
     [i1,i2,m]=regexp(fmt,'%[^FSL]*[FSL]','start','end','match'); args={};
     for l=1:length(i1)
        switch fmt(i2(l)) % S = module (subroutine)
           case 'F', fmt(i2(l))='s'; args{end+1}=stack.file;
           case 'S', fmt(i2(l))='s'; args{end+1}=stack.name;
           case 'L', fmt(i2(l))='d'; args{end+1}=stack.line;
           otherwise error('Wb:ERR','invalid fmt=`%s''',m{l});
        end
     end
     line=sprintf(fmt,args{:});
     return
  end

  if nargout==0, disp(stack)
  elseif nargout==1
   % default format: file>name:line
   % but drop name if equal to file (i.e. the actual function)
     if ~isequal(stack.file, [stack.name '.m']);
          line = sprintf('%s>%s:%d', stack.file, stack.name, stack.line);
     else line = sprintf('%s:%d', stack.file, stack.line); end
     if len>1, line=shortfstr(line,len); end
  else
     line=stack.line;
     name=stack.name;
  end

end

