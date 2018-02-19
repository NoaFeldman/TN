function i=isbatch
% Function: isbatch
%
%    Check wheter batchmode is set in get(0,'UserData').
%
% Wb,Sep06,07

  u=get(0,'UserData'); i=0;

% see also: (ismcc || isdeployed) for mcc-compiled routines!

  if isempty(u) || ~isfield(u,'batch')
   % consider jobs submitted to cluster as batch jobs
     if ~isempty(getenv('SGE_O_HOST')) ...
     && ~isempty(getenv('SGE_O_HOME')) i=1; end
     return
  end

  if ~isscalar(u.batch)
     wblog('ERR','Invalid batchflag in UserData !??');
     disp(u); return
  end

  i=u.batch;

end

