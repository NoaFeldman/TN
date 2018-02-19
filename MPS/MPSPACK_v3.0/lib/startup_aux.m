function startup_aux(pstr, istr)
% function startup_aux(pstr, istr)
% called upon startup through ~/bin/ml script.
% Wb,Feb01,06

% NB! by default, window title is set through ~/.bashrc :
% export PROMPT_COMMAND='echo -ne "\033]0;$(pwdprompt --bash-title)\007"'
% nevertheless, setenv PROMPT_COMMAND from within matlab has no effect!
% ! echo -ne "\033]0;$(pwdprompt --bash-title)\007"
% ! echo -ne "\033]0; GOTCHA :) \007"  # works like a charm!
% ==> see libx/cd.m which does exactly that! not that this is the very
%     1st command called in ~/bin/ml: $ML -r "cd $DIR; startup_aux()"
% Wb,Aug10,12

  I=mlinfo;

  nr={I.omp.num_threads, I.mkl.num_threads, I.cgi.num_threads};
  if any([nr{:}])
   % NB! matlab/2016a *unsets* MKL_NUM_THREADS at startup (!)
   % see ~/.matlab_setup.pl => MathWork communication
   % on MKL_NUM_THREADS, etc.! // Wb,Aug02,16
     if nr{1}>1 && (isempty(nr{2}) || nr{1}==nr{2})
        s=sprintf('%g! %g',nr{1},nr{3});
      % see MEX/tst_numthreads.m // Wb,Mar30,16
      % warning off MATLAB:maxNumCompThreads:Deprecated
      % n=maxNumCompThreads(NTH);
      % warning on MATLAB:maxNumCompThreads:Deprecated
        feature('numthreads',nr{1});
     else
        s=sprintf(',%g',nr{:}); s=s(2:end);
     end
     s=sprintf(', nthreads=[%s]',s);
  else s=''; end

  if I.cgi.verbose
  s=[s sprintf(', cgverbose=%g',I.cgi.verbose)]; end

  s=sprintf('pid=%g, %s%s', I.pid, get(0,'DefaultFigureRenderer'),s);

  if nargin==2
     fprintf(1,'\n   Matlab started by %s in %s\n   having %s\n\n',...
     repHome(pstr),repHome(regexprep(istr,'\<cd\> ','')),s);
  else t='';
     if I.status.ismcc,      t=[t ', mcc'];       end
     if I.status.isdeployed, t=[t ', deployed'];  end
     if I.status.isbatch,    t=[t ', batch'];     end
     if usejava('desktop'),  t=[t ', desktop'];   end
     if ~I.display,          t=[t ', nodisplay']; end
     if isempty(t), t='default'; else t=t(3:end); end

   % usejava('jvm')    %  1 even if matlab -nodisp, 0 if matlab -noj
   % usejava('awt')    %  0 if matlab -nodisp
   % usejava('swing')  %  0 if matlab -nodisp
   % usejava('desktop') : 0 if matlab is run in command-line mode(!)

     fprintf(1,'\n>> matlab environment: %s\n>> %s :: %s \n\n',...
     t,repHome(pwd),s);
  end

end

