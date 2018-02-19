function A=getvac(A,varargin)
% function A=getvac(A [,k])
%
%    reduce QSpace A to rank-2 Id tensor which only has
%    vacuum state / scalar representation with each index.
%    if k is specified, itag{k} will be used.
%
% Wb,Oct18,14

% outsourced to QSpace/getId.m

  A=getId(A,[],varargin{:});

end

