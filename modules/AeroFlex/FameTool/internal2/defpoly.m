function [x,y] = defpoly(i,gridHalf,gridPts,outline,bounds)
% Define polygon that will be used to find points inside the polygon

% Inboard line definition
pt1      = gridHalf(i,1:2);
pt2      = gridPts(i).coord(1:2);
[x,y]    = defline(pt1,pt2,1,bounds);
[xi, yi] = polypts(x, y, outline(:,1), outline(:,2));

% Outboard line definition
pt1      = pt2;
pt2      = gridHalf(i + 1,1:2);
[x,y]    = defline(pt1,pt2,2,bounds);
[xo, yo] = polypts(x, y, outline(:,1), outline(:,2));

% calculate intersection points with outline
pt1 = gridHalf(i,1:2);
pt2 = gridHalf(i + 1,1:2);

x = [xi;xo;pt1(1);pt2(1)];
y = [yi;yo;pt1(2);pt2(2)];

dx = gridPts(i).coord(1);
dy = gridPts(i).coord(2);

xd = x - dx;
yd = y - dy;

[t,~] = cart2pol(xd,yd);

[~,idx] = sort(t);

x = x(idx);
y = y(idx);

x = cat(1,x,x(1));
y = cat(1,y,y(1));
end