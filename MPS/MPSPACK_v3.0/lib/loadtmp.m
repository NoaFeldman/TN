function S=loadtmp(varargin)
% Function: S=loadtmp([opts,][usual arguments for matlab load])
%
%    load variables from temporary file
%    as data exchange from other MatLab session.
%
%    NB! if no return argument is specified,
%    the variables are loaded into the workspace.
%
% Options
%
%   '-<#>'    use temporary file with id <#> (#=1..9)
%
% See also save2tmp.m
% Wb,Nov12,09

  f=[ getenv('HOME') '/Matlab/tmp.mat' ];

  if nargin && ischar(varargin{1}) && ~isempty(regexp(varargin{1},'^-[0-9]+$'))
     tid=varargin{1}(2:end); varargin=varargin(2:end);

     q=str2num(tid); if q<1 || q>9
       error('Wb:ERR','\n   ERR invalid tid=%s',tid); end
     f=strrep(f,'tmp.mat',sprintf('tmp%g.mat',q));
  end

  if nargout==0
     cmd=['load ' f ' ' sprintf(' %s', varargin{:}) ];
     evalin('caller', cmd); 
  else
     S=load(f,varargin{:});
  end

end

