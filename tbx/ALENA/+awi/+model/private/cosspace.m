function x = cosspace(varargin)
%cosspace Generates a vector of points using cosine spacing.
%
% Inputs:
%   - Optional
%       + 'start'   : Start point of 'x'. 
%       + 'end'     : End point of 'x'.
%       + 'nPoints' : Number of points in 'x'.

%Parse inputs
p = inputParser;
addOptional(p, 'point1', 0, @(x)validateattributes(x, ...
    {'numeric'}, {'scalar', 'nonempty', 'finite', 'real', ...
    'nonnan'}, 'awi.model.CrossSection.cosspace', 'point1'));
addOptional(p, 'pointN', 1, @(x)validateattributes(x, ...
    {'numeric'}, {'scalar', 'nonempty', 'finite', 'real', ...
    'nonnan'}, 'awi.model.CrossSection.cosspace', 'pointN'));
addOptional(p, 'nPoints', 100, @(x)validateattributes(x, ...
    {'numeric'}, {'scalar', 'nonempty', 'finite', 'real', ...
    'nonnan'}, 'awi.model.CrossSection.cosspace', 'nPoints'));
addParameter(p, 'bRow', false, @(x)validateattributes(x, ...
    {'logical'}, {'scalar'}, 'awi.model.CrossSection.cosspace', 'bRow'));
parse(p, varargin{:});

%Grab values from parser
point1  = p.Results.point1;
pointN  = p.Results.pointN;
nPoints = p.Results.nPoints;

%Sort values in ascending order
if pointN < point1
    point1 = p.Results.pointN;
    pointN = p.Results.point1;
end

%Preallocate
x(nPoints, 1) = pointN;
x(1)          = point1;
midPoint      = (pointN - point1) / 2;

%Angle increment
dAngle = pi/(nPoints - 1);
%Vector of angles
angle = (dAngle : dAngle : (pi - dAngle));
%Array of x-values
x((2 : nPoints - 1)) = midPoint .* (1-cos(angle));

if p.Results.bRow
    x = x';
end

end