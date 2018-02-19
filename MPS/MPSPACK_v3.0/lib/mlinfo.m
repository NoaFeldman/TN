function I=mlinfo(varargin)
% function I=mlinfo()
% See also startup_aux.m
% Wb,Sep21,15

  getopt('init',varargin);
     dflag=getopt('--start');
  getopt('check_error');

  if dflag
     I.started=datestr(now);
     I.finished='';
  else
     I.date=datestr(now);
  end

  I.host={ hostname, hostid };
  I.pid=getpid();
  I.pwd=pwd;
  I.cputime=cputime;

  I.display=(usejava('awt') && usejava('jvm'));
  I.status=struct( ...
     'ismcc',     ismcc, ...
     'isbatch',   isbatch, ...
     'isdeployed',isdeployed, ...
     'isjava',  [ usejava('awt'), usejava('jvm') ],...
     'desktop',   usejava('desktop'), ...
     'Display', {{getenv('DISPLAY'), get(groot,'DefaultFigureRenderer')}} ...
  );

% usejava('jvm')    %  1 even if matlab -nodisp, 0 if matlab -noj
% usejava('awt')    %  0 if matlab -nodisp
% usejava('swing')  %  0 if matlab -nodisp
% usejava('desktop') : 0 if matlab is run in command-line mode(!)

% [i,s]=system('omp_threads.pl -x'); s=str2num(s);
% I.ncores=s(end);
% I.nslots=str2num(getenv('NSLOTS'));

% NB! maxNumCompThreads is deprecated with matlab >= R2009
% n=str2num(n); i=maxNumCompThreads(n);
% setuser(groot,'maxNumCompThreads',i);
% Wb,Feb10,11
  I.cgi=struct(... % CG_ info parameters
    'num_threads',str2num(getenv('CG_NUM_THREADS')), ...
    'verbose',    str2num(getenv('CG_VERBOSE')), ...
    'store',      getenv('RC_STORE')  ...
  );

% NB! matlab understands strings 'true/false' as 1/0 values
% => correctly treated in if-statements or also str2num()!
% Wb,Feb12,16
  I.omp=struct(...
    'num_threads',str2num(getenv('OMP_NUM_THREADS')), ...
    'max_active_levels', str2num(getenv('OMP_MAX_ACTIVE_LEVELS')), ...
    'dynamic',    getenv('OMP_DYNAMIC'), ...
    'nested',     getenv('OMP_NESTED'),...
    'schedule',   getenv('OMP_SCHEDULE'), ...
    'stack_size', str2num(getenv('OMP_THREAD_STACK_SIZE')) ...
  );

% https://software.intel.com/en-us/articles/recommended-settings-for-calling-intel-mkl-routines-from-multi-threaded-applications
% Wb,Feb12,16
  I.mkl=struct(...
    'num_threads',str2num(getenv('MKL_NUM_THREADS')), ...
    'num_threadD',str2num(getenv('MKL_DOMAIN_NUM_THREADS')), ...
    'dynamic',    getenv('MKL_DYNAMIC'), ...
    'nested',     getenv('MKL_NESTED'), ...
    'threading_layer',getenv('MKL_THREADING_LAYER') ...
  );

  I.jobid={ getenv('JOB_ID'), getenv('SGE_TASK_ID') };
  if isempty(I.jobid{2})
     if isempty(I.jobid{1}), I.jobid={};
     else I.jobid=I.jobid{1}; end
  end

end

