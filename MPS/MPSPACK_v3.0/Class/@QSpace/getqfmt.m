function f=getqfmt(A,varargin)
% function f=getqfmt(A [,opts])
%
%    get format string for a single Q{i}(j,:)
%
% Wb,Dec14,15

% see also display -> get_q_fmt(qtype,r,m)

  if isempty(A), f=''; return; end
  getopt('init',varargin);
     sep=getopt('sep',' ');
  getopt('check_error');

  d=0; f=getsym(A,'-c'); % if ~iscell(f), f={f}; end

  Q=abs(cat(1,A.Q{:}));

  for i=1:numel(f), s=f{i};
     if ~ischar(s), error('Wb:ERR',...
        '\n   ERR got invalid symmetry (string required) !?'); end
     if regexp(s,'^SU\d+$')
        r=str2num(s(3:end))-1; d=d+r; % r=N-1;
        f{i}=repmat('%X',1,r);
     elseif regexp(s,'^Sp\d+$')
        r=str2num(s(3:end))/2; d=d+r; % r=(2n)/2;
        f{i}=repmat('%X',1,r);
     else
        d=d+1;
        if any(Q(:,d)>10), f{i}='%3g'; else f{i}='%2g'; end
     end
  end

  if d~=size(Q,2), error('Wb:ERR',...
    '\n   ERR qset mismatch (len=%g/%g !?)',d,size(Q,2)); end

  f(2,1:end-1)={sep};
  f=[f{:}];

end

