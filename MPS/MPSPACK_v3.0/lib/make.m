function e=make(varargin)
% Function: make
%
%    compile source files given on input or default cc files instead.
%    Routine automatically moves to the MEX directory, but returns to pwd.
%
% Options
%
%   -f  force compilation
%   -u  udpate compilation if older than 1 day
%   -m  unix make using Makefile.mex
%
% See also mexall to recompile all existing mex files.
% Wb,Mar21,06

  if nargin==1 && isequal(varargin{1},'-e')
  vi make; return; end

  getopt('INIT',varargin);
     udays=getopt('-u',1);  % update within udays (default: 1)
     fflag=getopt({'-f','-B'});    % force compilation
     gflag=getopt('-g');    % debug option
     mflag=getopt('-m');
     tflag=getopt('-t');

     if getopt('-v'), vflag=2; % verbose
     elseif getopt('-V'), vflag=3; else vflag=0; end

     nolist=getopt('nolist');    % skip initial list
  varargin=getopt('get_remaining'); narg=length(varargin);

  mark=zeros(1,narg);
  for i=1:narg
     if ischar(varargin{i}) && ~isempty(varargin{i}) && varargin{i}(1)=='-'
     wblog('ERR','unrecognized option %s',varargin{i}); mark(i)=1; end
  end
  i=find(mark);
  if ~isempty(i), varargin(find(mark))=[]; narg=length(varargin); end

  if fflag, udays=-1; end

  if length(varargin)==0 || ~ischar(varargin{1})
  eval(['help ' mfilename]); return, end

% fls=dir('*.mexglx');
% fls=dir('*.mexa64');  % e.g. on 64bit processors such as duplo.

  R=getenv('MEXROOT');
  if isempty(R) || ~exist(R,'dir'), error('Wb:ERR',...
  'missing or invalid environmental variable MEXROOT'); end

  MEXD=[R,''    ];
  ODIR=[R,'/bin'];
  LDIR=[R,'/lib'];

  ext=mexext; tag=ext(end-2:end);
  
% COPTS: initial and trailing part of arguments handed over to mex
  COPTS={
     {'-outdir',ODIR}
     {'-lmwblas','-lmwlapack'}
     {'-lmpfr','-lgmp' }
  };
% COPTS={{'-outdir',ODIR}, {} };

% mex options
%  -g       debuggable
%  -inline  not possible (msg: Inlining not possible without C source files)

  LARP=['-L' R '/lib -larpack_' tag];

% NB! call -larpack BEFORE -llapack ( ARPACK requires LAPACK 2.0 !!)
%
%    otherwise, linking -llapack before -larpack takes lapack routines 
%    from -llapack (i.e. the CURRENT version, now lapack 3.0) and
%    IGNORES lapack routines also compiled into -larpack
%
% See also ARmake.inc in $ARPACK
% Wb,Aug27,08

  XFILES = {
   % NB! -lg2c is deprecated
   % arpack now recompiled with gfortran instead of f77
   % Wb,Feb09,10
     'eigs2Site', {'-lmwarpack'} % {LARP, '-lgfortran'} % '-lg2c'
     'wbeigs',    {'-lmwarpack'} % {LARP, '-lgfortran'} % '-lg2c'
  }; nx=size(XFILES,1);

  IMEX=['-I' MEXD ];

  fls_eigs={
    'eigs2SiteQS.cc'
  };

% compiles using make [-B] tst_omp :) # NB! do not use -m flag
% Wb,Jan25,11
% ********************************************************** %
% Wb,Jan30,12: -fopenmp included into mexopts.sh
% ********************************************************** %
% fls_tst_omp={
%  % NB! gcc-4.1 does not understand -fopenmp (only starting with gcc.4.2)
%  % NB! CFLAGS is for C-files! => use CXXFLAGS, instead!
%  % CLIBS CC=>CXXLIBS : libraries to link
%  % CC => CXX         : compiler to use
%  { 'tst_omp.cc','CC=gcc-4.3','CXXFLAGS=\$CXXFLAGS -fopenmp',...
%    'LDFLAGS=\$LDFLAGS -fopenmp' }
%  %  -Wall -x c++;  NB! -fopenmp => -pthread is implied
%  % tags: omp_ -lpthread
% };

% NB! use -t flag to show actual libraries linked
  fls_tst_gmpxx={ { 'tst_gmpxx.cc', '-lgmpxx -lgmp' }};
  fls_tst_mpfr_0={ { 'tst_mpfr_0.cc', '-lmpfr -lgmp' }};
  fls_tst_mpfr={ { 'tst_mpfr.cc','-lmpfr -lgmp' } };

  fls_mpfr2dec={ { 'mpfr2dec.cc','-lmpfr -lgmp' } };

% /amnt/dist64/DIR/matlab-R2011a/bin/glnxa64/libhdf5.so.6 /usr/lib/libhdf5.so
  fls_tst_ldbl={ { 'tst_ldbl.cc','-lmpfr -lgmp' } };

% GNU MP 4.1.2
% programs should be linked with "-lgmpxx -lgmp" libraries
  fls_tstmex={ { 'tstmex.cc', '-lmpfr -lgmp' }};

  fls_bicg={
    'bicg2SiteQS'
  };
  fls_bcg={
    'mpsOrthoQS'
    'bicg2SiteQS'
    'contractQS'
  };

  fls_tdm={
    'tdmNRG'
    'getSmoothTDM'
  };

  fls_fdm={
   {'fdmNRG_QS', '-lmpfr -lgmp'}
   %'getSmoothSpec'
  };

  fls_fgr={
    'fgrNRG'
   %'getSmoothSpec'
  };

  fls_nrg={
   {'NRGWilsonQS', '-lmpfr -lgmp' }
  % 'NRGWilsonQS'
  };

  fls_NRGWilsonQS={ {'NRGWilsonQS', '-lmpfr -lgmp' } };
  fls_fdmnrg={
   %'NRGWilsonQS'
   %'dmNRG'
   {'NRGWilsonQS', '-lmpfr -lgmp' }
   {'fdmNRG_QS', '-lmpfr -lgmp' }
   %'getSmoothSpec'
  };

  fls_getIdentityQS={ {'getIdentityQS', '-lmpfr -lgmp' } };
  fls_getid={
   {'getIdentityQS', '-lmpfr -lgmp' } % Wb,Nov21,12
  };

  fls_tst={ % 'phonebook.c', 'yprime.c', 'yprime2.c'
    'tstVEC.cc'
    'tstDATA.cc'
    'tstDGEMMc.cc'
    'tstmex.cc'
    'tst_mfp.cc'
  };

  fls_lib={
    'cell2matc.cc'
    'mat2cellc.cc'
    'printfc.cc'
    'getpid.cc'
    'dealstruct.cc'
    'uniquerows.cc'
  };

  fls_util = {
    'sparse2wb.cc'
    'kramers.cc'
    'getSmoothSpec.cc'
    'wbeig.cc'
    'wbrat.cc'
    'wbhist.cc'
    'getchar.cc'
  };

  fls_mps = {
    'bicgMatC'
    'contractIndex'
    'ivectorizeData'
    'matchIndex' % 'setup2SiteQ' 'svd2SiteQ'
    'vectorizeData'
    'getDimQS'
    'getDRangeQS'
    'permuteQS'
    'maxDiffQS'
    'plusQS'
    'getQDimQS'
    'getDimQS'
    'skipZerosQS'
    'permuteQS'
    'contractQS'
    'traceQS' % 'mpsTransferSpace'
    'mpsTensorProdQS'
    'mpsIsHConj'
    'mpsFullQS'
    'mpsFull2QS'
    'isDiagQS'
    'isIdentityQS'
    'isIdentityCG'
    'eigQS'
  };

  fls_cgs = {
    'compactQS'
    'contractQS'
    'getIdentityQS'
    'getSymStates'
    'compactQS'
    'contractQS'
    'isIdentityQS'
    'isIdentityCG'
  %{'NRGWilsonQS', '-lmpfr -lgmp' }
  %{'fdmNRG_QS', '-lmpfr -lgmp'}
  };

% fls_all = cat(1, fls_lib, fls_mps, fls_nrg, fls_bcg, fls_dma ); % fls_tst
% use 'mexall' instead

  fls_cgall = {
     'getSymStates'; 'compactQS'; 'contractQS'; 'getIdentityQS'
    {'NRGWilsonQS', '-lmpfr -lgmp' }
    {'fdmNRG_QS', '-lmpfr -lgmp' }
     'getCG'; 'wbtrace'
     'getDimQS'
     'getQDimQS'
     'getDRangeQS'
     'makeUniqueQS'
     'maxDiffQS'
     'normQS'
     'permuteQS'
     'plusQS'
     'skipZerosQS'
     'eigQS'
     'mpsFull2QS'
     'isDiagQS'
     'isIdentityQS'
     'mpsIsHConj'
     'isIdentityCG'
     'mpsOrthoQS'
     'timesElQS'
     'mpsTensorProdQS'
  };

  n=numel(varargin);
  if n && isequal(varargin{1},'all')
     wblog('<i>','\NUse ''mexall'' instead\N');
     return
  end

  for i=1:n
     f1=['fls_' varargin{i}];
     if exist(f1,'var')
        eval(['fls=' f1 ';']); varargin{i}=fls;
        if vflag
           if ~iscell(varargin{i})
                fprintf(1,'\n   %s => %s:\n',varargin{i},f1);
         % else disp(strhcat(varargin{i});
           end
         % disp(fls);
        end
     else varargin{i}={varargin{i}};
     end
  end
% varargin
  fls=cat(1,varargin{:}); n=numel(fls);

% NB! can sort cell array of strings
% but must ensure, there is no further cell substructure.
% Wb,Aug13,12
  if n>1
     for i=1:n, if iscell(fls{i}), i=i-1; break; end, end
     if i==n, fls=unique(fls); end
  end

% keyboard; return

  nf=length(fls); if ~vflag && nf>1, vflag=1; end
  if vflag>2, COPTS{1}={'-v',COPTS{1}{:}}; end

  cdir=pwd; cd(MEXD); mcount=0; e=zeros(1,nf);

  if vflag && ~nolist
     m=0; ff=fls; % inl(1);
     for i=1:numel(ff), if iscell(ff{i}), ff{i}=strhcat(ff{i}); m=m+1; end, end
     if m, fprintf(1,'   %s\n',ff{:}); 
     else fprintf(1,'%s ',[char(10) '=>'],fls{:},char(10)); end
     inl(1);
  end

  for k=1:nf, f=fls{k}; dl=0;

    % NB! -largeArrayDims required by mxGet(Ir|Jc)
    % disabled MX_COMPAT_32 in wblib.h (no longer applicable for 64bit systems)
    % Wb,Jan07,11
    % copts={}; % o={'-compatibleArrayDims'};
    % -silent for extra messages of matlab 2014b // Wb,Jan24,15
      copts={'-largeArrayDims'}; % ,'-silent'

    % mexlib.h ExitMsg() is sensitive to #define DBSTOP
    % => make compiler stop at error messages
      if gflag
         copts={'-g',copts{:}};
         if ~isempty(getenv('DBSTOP')), copts{end+1}='-DDBSTOP'; end
       % kill standard MatLab session otherwise!
      end

      if iscell(f)
         if vflag
            fprintf(1,'=> %s%s\n\n',f{1},sprintf(' %s',f{2:end}));
         end
         copts={copts{:},f{2:end}}; f=f{1}; dl=1;
      else for j=1:nx
         if ~isempty(findstr(f,XFILES{j}))
         copts={XFILES{j,2}{:},copts{:}}; % keep this order! tags: BEFORE
         end
      end,end

      copts={COPTS{1}{:}, copts{:}, COPTS{2}{:}};

      if isempty(findstr(f,'.c'))
       % skip extension (if any) and try to find .c or .cc file on the path
         i=find(f=='.'); if ~isempty(i), f=f(1:i(end)-1); end
         full=which([f '.cc'],'-all'); if isempty(full),
         full=which([f '.c'],'-all'); end
      else
         full=which(f,'-all');
      end

      n=length(full);
      if n==1, full=full{1};
      else
         if n>1, s='several files found';
         else s='unable to find source file'; end
         s=sprintf([s ' for >%s<'],fls{k});

       % error('Wb:ERR',s);
         wblog('ERR',s); % keyboard
         if n, disp(strvcat(full{:})); end
         e(k)=1; continue
      end

      q=regexp(copts,'mpfr'); if isempty([q{:}])
         if str2num(evalc(['! egrep -c CGC_QSPACE "'  full '"']))
          % fprintf(1,'   ... adding %s\n',strhcat(COPTS{3}));
            copts=[copts, COPTS{3}];
         end
      end

      [p,f,x]=fileparts(full); fext=[f '.' ext];
      fulx=which(fext,'-all');

      n=length(fulx);
      if n==1
         fulx=fulx{1};
         pout=fileparts(fulx); % path to mex file

         if ~isequal(pout,ODIR) && ~isequal(pout,LDIR) && ~isequal(pout,MEXD)
         wblog('NB!','mex file `%s''\nwritten to `%s''',fext,pout); end

         setopts(copts,{'-outdir',pout},'-q');

      elseif n<1
         fulx=''; pout='';

       % c-files in LDIR=MEX/lib are also compiled into that directory
         if isequal(p,MEXD), s=ODIR;
         else s=p; end % overwrite existing file! (including LDIR)

         setopts(copts,{'-outdir',s},'-q');
      else
         s=sprintf('several mex-files found for >%s<',fls{k});

       % error('Wb:ERR',s);
         wblog('ERR',s); % keyboard
         disp(strvcat(fulx{:})); e(k)=1; continue
      end

      if ~isempty(p) && ~isequal(p,MEXD)
        % current directory is MEXD
          copts={copts{:}, IMEX};
          cd(p)
      end

    % fprintf(1,'%30s %-16s %-16s %d %d',' ', fext, full, ...
    % ftime(fext) < ftime(full), exist(fext,'file'))

      if fflag || isempty(pout)  ...
         || ftime(fulx) < ftime(full) || ftime(fulx) < now-udays

         if vflag || dl
            s=sprintf('-> %2d/%d: mex [%s] %s',k,length(fls),tag,full);

            if ~isempty(copts)
            s=[s, sprintf(' %s', copts{:})]; end
            s=regexprep(s,'-lmw\w+\s*',''); % skip lapack blas entries
            s=strrep(s,'-outdir','=>');

            s=strrep(strrep(s,'/amnt',''),MEXD,'$MEX');
            s=regexprep(s,'\s*-largeA\w*','');

            fprintf(1,'%s ...',s);
         end

         mcount=mcount+1; s=sprintf(' %s',copts{:});
         if ~isempty(findstr(s,'-L'))
            s=regexprep(s,'(-out|-L)',[char(10) '   $1']);
            fprintf(1,'\nmex %s\n\n',...
            regexprep([full ' [' tag ']',s],'/[^ ]*MEX','$MEX'));
         end

       % setopts(copts,{'-outdir'}); copts=copts([1 3]); % tags: COPTS
       % strvcat(full, copts{:})

         if vflag>2
            fprintf(1,'\n\n  mex %s\n\n',strhcat(copts{:},'-s',char(10)));
         end

         if mflag, fprintf('\n '); od=''; oc=''; ll={};
            for i=1:length(copts)
                if isequal(copts{i},'-outdir'), od=copts{i+1};
                elseif ~isempty(regexp(copts{i},'^-[ILl]'))
                   ll{end+1}=copts{i}; % -L.. -I.. -lmwblas ...
                end
            end

            if fflag, oc='-B'; end
            cmd=[ 'export MXC_ODIR=''' od '''' 10 ...
                  'export MXC_LIB=''' strhcat(ll) '''' 10 ...
                  'export LOG_FLAG=1' 10 ...
            'make ' oc ' -f Makefile.mex ' strrep(full,'.cc','')];

            if tflag, fprintf(1,'\ncmd: %s\n\n',cmd); else
               i=system(cmd); if i, error('Wb:ERR',''); end
            end

         else
            try
             % temporary fix due to system upgrade // Wb,Jan09,15
             % copts{end+1}='-I/usr/include/x86_64-linux-gnu';

               if tflag
                  fprintf(1,'\nmex(%s,%s)\n\n',full,strhcat(copts)); 
               else
                  if vflag
                     i=find(findstrc(copts,'outdir'));
                     if ~isempty(i)
                        fprintf(1,'\n    ... %s : %s',copts{i},copts{i+1});
                        i=i+2:numel(copts);
                     else
                        i=1:numel(copts);
                     end
                     fprintf(1,'\n    ... %s\n',strhcat(copts(i))); 
                  end
                  mex(full,copts{:}); 
               end
            catch
               fprintf(1,'\n%s\n%s\n\n',full, sprintf('%s ',copts{:}))
             % full, strvcat(copts{:})
               cd(cdir)
               if nargout, e(k)=1; else 
               rethrow(lasterror); % yeah right B) ...
               end
            end
         end
         fprintf(1,'\n');
      else
         if vflag || length(fls)<=1
            fprintf(1,'   %s is up to date.\n', ...
            strrep(strrep(full,'/amnt',''),MEXD,'$MEX'));
         end
      end
  end

  if vflag
     if mcount==0 && sum(e)==0
     fprintf(1,'\n   All %d file(s) uptodate.\n',nf);
     elseif length(fls)>0, inl 1; end
  end

  cd(cdir)
  if ~nargout, clear e; end

end

