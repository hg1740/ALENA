function Euler = Quat2Euler(quat)

phi   = atan2(2*(-quat(1)*quat(2) + quat(3)*quat(4)),1-2*(quat(2)^2+quat(3)^2));
theta = asin( 2*(-quat(1)*quat(3) - quat(4)*quat(2)));
psi   = atan2(2*(-quat(1)*quat(4) + quat(2)*quat(3)),1-2*(quat(3)^2+quat(4)^2));

Euler = [phi;theta;psi];