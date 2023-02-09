function xskew = skewmat(x)
%skewmat Matrix version of 'skew'.
%
% Returns a 3D skew matrix of a [3, n] matrix.
%
% See also: skew

%Preallocate
xskew = zeros(3, 3, size(x, 2));

%Assign data
xskew(1, 2, :) = -x(3, :);
xskew(1, 3, :) =  x(2, :);
xskew(2, 1, :) =  x(3, :);
xskew(2, 3, :) = -x(1, :);
xskew(3, 1, :) = -x(2, :);
xskew(3, 2, :) =  x(1, :);

end