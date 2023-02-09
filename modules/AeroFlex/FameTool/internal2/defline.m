function [x,y] = defline(pt1,pt2,refpoint,bounds)

x1 = pt1(1);
y1 = pt1(2);
x2 = pt2(1);
y2 = pt2(2);

y = bounds(1):1:bounds(2);

m = (y2 - y1) / (x2 - x1);
m =  - 1 / m;

if m ==  0 && refpoint ==  1
    x = 0:1:50;
    y = repmat(y1,size(x));
elseif m ==  0 && refpoint ==  2
    x = 0:1:50;
    y = repmat(y2,size(x));
elseif refpoint ==  1
    c = y1 - m * x1;
    x = (y - c) / m;
else
    c = y2 - m * x2;
    x = (y - c) / m;
end




end
