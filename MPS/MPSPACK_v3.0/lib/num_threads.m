function n=num_threads(n)
% function n=num_threads(n)
%
%    (Re)set number of maximum active computational threads.
%
% Wb,Feb21,17

% undocumented matlab (yet recommended by matlab customer support)
  n0=feature('numthreads');
  if ~nargin; n=n0; return; end
% if n==n0, return; end

  feature('numthreads',n);
  
  s=sprintf('%g',n);
  setenv('OMP_NUM_THREADS',s);
  setenv('MKL_NUM_THREADS',s);

  setenv('MKL_DYNAMIC','false');
  setenv('OMP_NESTED','true');

% it appears MKL_NUM_THREADS is *unset* when matlab starts
% irrespective to what is set above, whereas MKL_DOMAIN_NUM_THREADS
% is newly introduced (see ~/.matlab_setup.pl // Wb,Feb08,17
% setenv('MKL_DOMAIN_NUM_THREADS',s);

% setenv('QSP_NUM_THREADS',s);

end

