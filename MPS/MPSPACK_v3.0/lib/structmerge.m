function ss=structmerge(varargin)
% Function: ss=structmerge(struct1, struct2, ...)
%
%    merge given set of structures
%    * structures can be arrays with equal dimensions each
%    * cell arrays of structures are automatically converted to structure arrays
%
% Wb,May16,07

  if nargin<2
     eval(['help ' mfilename]);
     if nargin, error('Wb:ERR','invalid usage'); end
  end

  s1=size(varargin{1});
  m=numel(varargin{1});

  for k=1:nargin
     if iscell(varargin{k})
     varargin{k}=reshape(cat(1,varargin{k}{:}),s1); end

     if ~isstruct(varargin{k})
        eval(['help ' mfilename]);
        error('Wb:ERR','all input arguments must be structures!');
     end
     if ~isequal(size(varargin{k}),s1)
        eval(['help ' mfilename]);
        error('Wb:ERR','size mismatch of input arguments');
     end
  end

% make sure fields are unique
  f=cell(1,nargin); nf=zeros(1,nargin);
  for k=1:nargin, f{k}=fieldnames(varargin{k}); nf(k)=length(f{k}); end
  ff=cat(1,f{:})';

  [x,I,d]=uniquerows(strvcat(ff{:}));
  if any(d>1)
   % wblog(1,'WRN','overwriting existing fields by later entries (%g)',sum(d-1));
   % NB! keep the field order according to which field is first appearing!
     for i=1:length(I), I{i}=I{i}(1); end
     ff=ff(sort(cat(2,I{:})));
  end
  
  ff(2,:)={[]};

  ss=struct(ff{:}); ss=ss(ones(s1));
  for i=1:m, S=ss(i);
      for k=1:nargin, s=varargin{k}(i);
         for j=1:nf(k), fn=f{k}{j};
         S=setfield(S,fn,getfield(s,fn)); end
      end
      ss(i)=S;
  end

end

