function [xCoords, zCoords] = scaleProfileData(xCoords, zCoords, xzData, dim)
%scaleProfileData Returns the coordinates to their original dimensions
%using the scale/translate data in 'xzData'.

nSection = size(xzData, 1);

if nargin < 4
   if size(xCoords, 1) == nSection
       dim = 2;
   elseif size(xCoords, 2) == nSection
       dim = 1;
   else 
       dim = [];
   end
end

if dim == 2
    n   = size(xCoords, 2);
    xSc = repmat(xzData(:, 3), [1, n]);
    zSc = repmat(xzData(:, 4), [1, n]);
    x0  = xzData(:, 1);
    z0  = xzData(:, 2);    
elseif dim == 1    
%     xCoords = xCoords';
%     zCoords = zCoords';    
    n   = size(xCoords, 1);
    xSc = repmat(xzData(:, 3)', [n, 1]);
    zSc = repmat(xzData(:, 4)', [n, 1]);
    x0  = xzData(:, 1)';
    z0  = xzData(:, 2)';
else 
    error('Unable to proceed. Unsure how to scale profiles.');
end

xCoords = xCoords .* xSc + x0;
zCoords = zCoords .* zSc + z0;

end