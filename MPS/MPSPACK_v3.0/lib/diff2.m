function [xd,yd]= diff2(varargin)
% Function: diff2 - modified diff() routine
%
%    returns the k-th derivative of arr columns given
%    where the first columns is treated as x vector;
%    an  array of dimension (m-k) x n gets returned!
%
% Usage
%
%    yd = diff2(yd)
%    [xd,yd] = diff2(xd,yd)
%
% Options
%
%    'k',..   k-th derivative (1)
%    'len'    preserve length
%
% Wb,Jun19,01

  if ~nargin, eval(['help ' mfilename]); return; end

  if nargin>1 && isnumeric(varargin{2})
     xd=varargin{1};
     yd=varargin{2};
     varargin=varargin(3:end);
  else
     xd=[];
     yd=varargin{1};
     varargin=varargin(2:end);
  end

  getopt ('init', varargin);
     k    =getopt('k',1);
     lflag=getopt('len'); % keep length of vector
  varargin=getopt('check_error');

  if ~isempty(xd)
     s=size(yd); n=numel(xd);
     if ~isvector(xd) || numel(s)>2 || all(s~=n)
        error('Wb:ERR', ...
       'size mismatch (%d; %s)', length(xd), sizestr(yd));
     end

   % take data along columns
     if size(xd,1)<size(xd,2), tx=1; xd=xd.'; else tx=0; end
     if s(1)~=n, ty=1; yd=yd.'; else ty=0; end
     m=size(yd,2);

     if ~isvector(xd) || length(xd)~=size(yd,1)
        error('Wb:ERR', ...
        sprintf('Size mismatch (%d/%d)', length(xd), length(yd)));
     end

     for i=1:k
        dx=1./diff(xd);
        yd=diff(yd) .* dx(:,ones(1,m));

        if lflag
         % calculate derivative on xi not midway between xi,i+1
           yd=[yd(1,:); 0.5*(yd(1:end-1,:)+yd(2:end,:)); yd(end,:)];
        else
         % set xi to midway between xi,i+1
           xd=0.5*(xd(1:end-1)+xd(2:end));
        end
     end

     if tx, xd=xd.'; end
     if ty, yd=yd.'; end

     if nargout==1, xd=yd; end
  else
   % take data along columns
     if size(yd,1)<size(yd,2), ty=1; yd=yd.'; else ty=0; end

     for i=1:k
        yd=diff(yd);

      % calculate derivative on xi not midway between xi,i+1
        if lflag
        yd=[yd(1,:); 0.5*(yd(1:end-1,:)+yd(2:end,:)); yd(end,:)]; end
     end

     if ty, yd=yd.'; end

     xd=yd; clear yd; % only one variable to return
  end

end

