% Load a NASTRAN file and create a beam struct
% Calculates rotation Matrix given Gibbs-Rodrigues parameters g
function R = Rmat(gin)
g = tan(gin/2);
%g = tan(gin/2);
%g=gin;
%R = zeros(3,3);
%crossmg=crossm(g);
%R = eye(3) + (4 / (4 + g'*g)) .* (crossmg + 0.5.* crossmg * crossmg);
%R = eye(3) + (2 / (1 + g'*g)) .* (crossmg + crossmg * crossmg);
R = [1 0 0 ;0 1 0;0 0 1] + (2 / (1 + g(1)^2+g(2)^2+g(3)^2)) * ([-g(3)^2-g(2)^2,g(2)*g(1)-g(3),g(3)*g(1)+g(2);g(1)*g(2)+g(3),-g(3)^2-g(1)^2,g(3)*g(2)-g(1);g(1)*g(3)-g(2),g(3)*g(2)+g(1),-g(2)^2-g(1)^2]);
