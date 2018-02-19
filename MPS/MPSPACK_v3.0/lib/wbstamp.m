function s=wbstamp(varargin)
% function s=wbstamp(varargin)
% Wb,Jun18,07

  getopt('INIT',varargin);
     lflag=getopt('-l'); if ~lflag && getopt('-L'), lflag=2; end
     fflag=getopt('-f');
     sflag=getopt('-s'); % short (without hour indicator)
     iflag=getopt('-init');
     cflag=getopt('-clear');
  getopt('check_error');

  if cflag, getuser(groot,'wbstamp','-rm'); end
  if ~iflag
     s=getuser(groot,'wbstamp'); % Wb,May15,13
     if ~isempty(s) && ischar(s) && ~fflag
        wblog(' * ','using wbstamp=''%s''',s);
        return
     end
  end

% if ~isempty(getenv('JOB_ID'))
%    s=strtrim(evalc('!wbstamp.pl')); % remove newlines
%    return
% end


  if sflag, h=''; else
   % add hour indicator through single extra char
     c=clock; h='a'+c(4);
  end

  s=['Wb', datestr(now,'yymmdd'), h];

% s=['Wb' datestr(now,'yymmdd_HHMMSS')]; interchanges HH with dd !_)$(!_
  if lflag
     if lflag==1,
          s=[s '_' datestr(now,'HHMM')];
     else s=[s '_' datestr(now,'HHMMSS')]; end
  end

  if iflag
     wblog(' * ','setting persistent wbstamp=''%s''',s);
     setuser(groot,'wbstamp',s);
     if ~nargout, clear s, end
  end

end

