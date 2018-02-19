function ss=catc(ic,varargin)
% function ss=catc(d,varargin)
%
%    similar to cat, but specifically meant for structures
%    adjust number of fields / field order if required.
%
% Wb,Jan28,10

  if ~nargin || ic~=1 && ic~=2
     eval(['help ' mfilename]);
     if nargin || nargout, error('Wb:ERR','invalid usage'), end, return
  end

  try
     ss=cat(ic,varargin{:});
  catch
   % structures not compatible
     n=numel(varargin); fn=cell(n,1); fs={}; m=0;
     for i=1:n, if isempty(varargin{i}), continue; end
        if ~isstruct(varargin{i})
        error('Wb:ERR','invalid usage (set of structures expected)'); end
        fn{i}=fieldnames(varargin{i});
      % get started with structure with most fields
        if length(fn{i})>m, fs=fn{i}; m=length(fs); end
     end
   % complete fieldnames fs
     for i=1:n, if isempty(fn{i}), continue; end
        m=length(fn{i});
        s=strvcat(fn{i}{:},fs{:});
        [x,j]=setdiff(s(1:m,:),s(m+1:end,:),'rows');
        if ~isempty(j), fs={fs{:},fn{i}{j}}; end
     end

     fs=reshape(fs,1,[]); fs(2,:)={[]};
     ss=struct(fs{:}); fs=fs(1,:); nf=numel(fs);

   % adjust fields of all structures
     for i=1:n, q=varargin{i};
        for j=1:nf, if ~isfield(q,fs{j})
           q=setfield(q,{1},fs{j},{1},{[]});
        end, end
        varargin{i}=orderfields(q,ss);
     end

     ss=cat(ic,varargin{:});
  end

end

