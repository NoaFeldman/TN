function tstfun(A)
% function tstfun(A)
% Wb,Feb23,11

  wblog('1'); whos A

% BEHAVIOUR WITHIN CLASS METHODS
% whenever a class method requires the functionality of the overloaded
% subsref or subsassign, it must call the overloaded methods with
% FUNCTION CALLS rather than using operators like '()', '{}' or '.'

  A(1) % call built-in MatLab subsref() method (!) -> return the whole A
  subsref(A, struct('type', '()', 'subs', {{1}}) )
% subsref(A, struct('type', '{}', 'subs', {{1}}) )
% A{1} % call built-in MatLab subsref() method (!) -> produces an error

  wblog('2');

  helpthis
  tstfun_aux(A)

end

function tstfun_aux(A)
% function tstfun_aux(A)
% Wb,Feb23,11

% access help for *this subroutine by typing
% >> help @QSpace/tstfun>tstfun_aux

  s=dbstack; s=s(1);
% f=str2func(s.name)
% help([ '@' class(A) '/' s.file '>' s.name ])
  helpthis

end

