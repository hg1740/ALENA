function w = getGustSensitivity(X,Y,Z,IN)

X_tmp = Y;
Y     = X;
X     = X_tmp;

H       = IN.H;
yg      = IN.yg;
wgmax   = IN.wgmax;

% if Y < -yg-X && Y > -yg-2*H-X
%     wx = 0;
%     wy = 0;
%     wz = wgmax*0.5*(1-cos((Y+yg+X)*pi/H));
% else
%     wx = 0;
%     wy = 0;
%     wz = 0;
% end

% if Y < -yg
%     wx = 0;
%     wy = 0;
%     wz = cos(X/5)*sin((Y+yg)/10)*wgmax;
% else
%     wx = 0;
%     wy = 0;
%     wz = 0;
% end

% if Y < -yg && Y > -yg-2*H && abs(X)<25
%     wx = 0;
%     wy = 0;
%     wz = wgmax*0.5*(1-cos((Y+yg)*pi/H))*X/16;%wgmax*0.5*(1-cos((Y+yg)*pi/H))*cos(X*pi/50);
% else
%     wx = 0;
%     wy = 0;
%     wz = 0;
% end

theta = IN.theta;

wx = 0;
wy = 0;
wz = wgmax*0.5*(1-cos((Y+yg)*pi/H));

wy = wz*sin(theta);
wz = wz*cos(theta);

w = [wx;wy;wz];