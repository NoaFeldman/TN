function [g0,a0,g1]=getConductance0(om,a0,varargin)
% function [g0,a0]=getConductance0(om,a0 [,T,gfac,...])
%
%    Get conductance by folding discrete a0 data with
%    derivative of Fermi distribution (-df/dw).
%
% Options
%
%   '-f'  got functional data (rather than discrete data set)
%
% outsourced from MEX//rnrg_gg.m ; tags: conductivity
% Wb,Jul08,13

  if nargin<2 || ~isvector(om)
     helpthis, if nargin || nargout
     error('Wb:ERR','invalid usage'), end, return
  end
  e=0;fflag=0;

  if isempty(varargin) || ischar(varargin{1})
     getbase('Idma');
     if isfield(Idma,'T'), T=Idma.T; else e=e+1; end
  else
     T=varargin{1}; varargin=varargin(2:end);
  end

  if isempty(varargin) || ischar(varargin{1})
     getbase('Inrg');
     if isfield(Inrg,'param') && isfield(Inrg.param,'Gamma')
        gfac=pi*Inrg.param.Gamma;
     else e=e+1; end
  else
     gfac=varargin{1}; varargin=varargin(2:end);
  end

  if e, error('Wb:ERR','\n   ERR invalid usage'); end

  if numel(varargin)
     getopt('init',varargin);
        fflag=getopt('-f');
     getopt('check_error');
  end

  if size(om,2)~=1, om=om'; end

  beta=1/T;
  fw=(beta/2)./(1+cosh(beta*om)); % derivative of Fermi-function

  if fflag
     dw=[ om(2)-om(1); 0.5*(om(3:end)-om(1:end-2)); om(end)-om(end-1) ];
     fw=fw.*dw;
  end

  g0=gfac*(fw'*a0);

  if nargout>1
     if nargout>2
      % for consistency with rnrg_gg.m
        i=find(abs(om)<=1);
        g1=gfac*(fw(i)'*a0(i,:));
     end

   % get A(omega=0) value
   % i=find(abs(om)<=dE); a0=sum(a0(i,:)/(2*dE));
   % i=find(abs(ox)<=dE); a0=mean(ax(i,:));
     dE=T/10; fw=exp(-(om/dE).^2) / (sqrt(pi)*dE); % Gaussian around om=0
     a0=fw'*a0;
  end

end

