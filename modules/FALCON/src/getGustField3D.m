function OUT = getGustField3D(IN)

% [x_coord,y_coord,z_coord] = meshgrid(-25: 1: 25,25:-1:-250,-25:1:25);
% [x_field,y_field,z_field] = deal(zeros(size(x_coord)));
% 
% H                  = 5;
% y_start            = 10;
% wg                 = 5;

[x_coord,y_coord,z_coord] = meshgrid(IN.x_range,IN.y_range,IN.z_range);
[x_field,y_field,z_field] = deal(zeros(size(x_coord)));

H       = IN.H;
y_start = IN.y_start;
wg      = IN.wg;

% x_field(y_coord>-(2*H+y_start) & y_coord<-y_start) = wg*0.5*(1-cos((y_coord(y_coord>-(2*H+y_start) & y_coord<-y_start)-y_start)*pi/H));
z_field(y_coord>-(2*H+y_start) & y_coord<-y_start) = wg*0.5*(1-cos((y_coord(y_coord>-(2*H+y_start) & y_coord<-y_start)-y_start)*pi/H));%.*x_coord(y_coord>-(2*H+y_start) & y_coord<-y_start)/25;

% x_range = -25:1:25;
% streamline(x_coord,y_coord,z_coord,x_field,y_field-25,z_field,x_range,25*ones(length(x_range),1),zeros(length(x_range),1))
% quiver3(x_coord,y_coord,z_coord,x_field,y_field-25,z_field)

OUT.x_coord = x_coord;
OUT.y_coord = y_coord;
OUT.z_coord = z_coord;
OUT.x       = x_field;
OUT.y       = y_field;
OUT.z       = z_field;