function clr(varargin)
% Function: clr [OPTS] - clear entire workspace
% Options:
%   -k  keep figures, classes, debug stops
% Wb,Jan05,07

  try
  getopt ('init', varargin);
     keep = getopt ('-k');
     figf = getopt ('-f');
     eflag= getopt ('-e');  % call subsequent 'dbstop if error'
     rflag= getopt ('-r');  % reset terminal // Wb,Aug27,15
  varargin=getopt('check_error');
  catch
     return
  end

% ----------------------------------------------------------------------- %
% Wb,Ap21,06 - clear classes when construtors is changed
% (rather than restarting matlab!)
% NB! clear functions also clears dbstop's

  if ~keep
     if ~figf
     evalin('caller','close all'); else untagf -a -q; end
     evalin('caller','clear functions; dbclear all; clear classes');
  end

  evalin('caller','clear all; clear global');

  if eflag, dbstop if error; end
  if rflag, system('reset'); end

  if keep || figf, untagf -a -q; return; end
  if rflag, return; end

% ----------------------------------------------------------------------- %
% tic

% if isequal(getenv('COLORTERM'),'gnome-terminal')
     s=repmat('%',1,80);
     fprintf(1,'\n%% %s %%\n%% %25s%-55s %%\n%% %s %%\n',...
     s,'','C  L  E  A  R      A  L  L',s); clc
     return
% end

% LINES is set via resize in ~/bin/ml script
  L=getenv('LINES'); if ~isempty(L), L=str2num(L); else L=30; end

% p={'',  -2};
% p={'clr_aux2.txt',  0};
  p={'clr_line.txt',  0};
% p={'clr_eagle.txt',14};
% p={'clr_aux.txt',  14};
% p={'',             -2};
% p={'clr_aux0.txt',  4};
% p={'clr_auxx.txt',  8};

  nl=10; n=ceil((L-p{2})/2)-4;

  disp([10 '%% CLEAR ALL ' repmat('%',1,73) 10 10]);

  disp(strvcat(' ', nl(ones(1,n))))

  if ~isempty(p{1}), type(p{1}); end

% disp(char('='+zeros(1,80)))
  disp(strvcat(' ', nl(ones(1,n))))

% pause(max([0.5-toc, 0]));

  clear; clc

end

