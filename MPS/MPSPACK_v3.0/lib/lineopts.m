function opts=lineopts(lh)
% Functions: opts=lineopts(lh)
% 
%    return lineopts from given line handle as cell array
%    if lh is an array, then only the line options of the
%    first element are returned.
% 
% Wb,Nov22,07 ; Wb,Aug04,16

  if nargin~=1 || ~ishandle(lh)
     helpthis, if nargin || nargout
     error('Wb:ERR','invalid usage'), end, return
  end

  if isempty(lh), error('Wb:ERR','got empty handle'); end
  if ~all(isline(lh(:))), get(lh,'Type')
     error('Wb:ERR','got non-line handlese');
  end

  lo={'Color','LineStyle','LineWidth',...
      'Marker','MarkerSize','MarkerEdgeColor','MarkerFaceColor',...
      'DisplayName','Tag' };

  lo(2,:)=get(lh(1),lo);
  opts=reshape(lo,1,[]);

end

