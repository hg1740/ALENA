function [xNorm, zNorm, xzData] = normaliseProfileData(x, z)
%normaliseCrossSection Normalises the profile coordinates between [0 : 1]
%and centres the cross-section on the (0.5, 0.5) coordinate.
%
% Returns the normalised coordinates and a matrix of data, each column
% represents the offset and scale data for the (x,z) axes.
%
% e.g. xzData = [xMin, zMin, xScale, zScale]

xMin   = min(x, [], 2);
xMax   = max(x, [], 2);
zMin   = min(z, [], 2);
zMax   = max(z, [], 2);

xScale = (xMax - xMin);
zScale = (zMax - zMin);

xNorm  = (x - xMin) ./ xScale;
zNorm  = (z - zMin) ./ zScale;

xzData = [xMin, zMin, xScale, zScale];
end