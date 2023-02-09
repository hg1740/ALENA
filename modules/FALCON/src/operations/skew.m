function xskew = skew(x)
%skew Returns the skew matrix of a vector. The skew matrix is analogous to
%the cross product of two [1,3] vectors. CHECK THIS DEFINITION!!!
%
% xskew = [ ...
%      0   , -x(3),  x(2); ...
%      x(3),  0   , -x(1); ...
%     -x(2),  x(1),   0  ];
xskew = [0, -x(3), x(2); x(3), 0, -x(1); -x(2), x(1), 0 ];

end
