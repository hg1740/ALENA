% Calculates rotation Matrix G given Gibbs-Rodrigues parameters g
function G = Gmat(gin)
g = tan(gin/2);
% G = zeros(3,3);
%G = (4 / (4 + g'*g)) .* (eye(3) + 0.5*crossm(g));

G = (1 / (1 + g(1)^2+g(2)^2+g(3)^2))*([1,-g(3),g(2);g(3),1,-g(1);-g(2),g(1),1]);
