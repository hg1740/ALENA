function [xU, xL, zU, zL, etaQ2D] = interpolateProfileData(x, z, eta, etaQ, n)
%interpolateProfileData Interpolates the normalised aerofoil coordinates
%defined at 'eta' to the new spanwise position 'etaQ'.
%
% Inputs:
%   - 'x'    : Profile x-coordinates at positions 'eta'.
%   - 'z'    : Profile z-coordinates at positions 'eta'.
%   - 'eta'  : Nondimensional spanwise positions of profiles.
%   - 'etaQ' : Nondimensional spanwise positions where the profiles will be
%              interpolated to.
%   - 'n'    : Number of points in the x-direction for the upper and lower
%              surfaces of the profile.

if nargin < 5
    n = size(x, 1) / 2;
end

%Split into upper and lower surfaces
%   - Prevents the "non-monotincally increasing" error when
%     interpolating
XUpper_ = x(:, 1 : n);
XLower_ = x(:, n + 1 : end);
ZUpper_ = z(:, 1 : n);
ZLower_ = z(:, n + 1 : end);

%Check the x-coords are defined at the same points
%   - This means we are using the same grid to perform the 2D interpolation
%   when we create the mesh grid using only the first set of x-coords
if or( ...
        any(any(abs(diff(XUpper_, [], 1)) > 1e-10)), ...
        any(any(abs(diff(XLower_, [], 1)) > 1e-10)))   
    error(['Unable to interpolate cross-sections that ', ...
        'are not defined for the same normalised '     , ...
        'x-coordinates.']);
end

%Make meshgrids for interpolation
[xU_old, e]  = meshgrid(XUpper_(1, :), eta');
[xL_old, ~]  = meshgrid(XLower_(1, :), eta');
[xU, etaQ2D] = meshgrid(XUpper_(1, :), etaQ);
[xL, ~]      = meshgrid(XUpper_(1, :), etaQ); 

%2D interpolation
zU = interp2(xU_old, e, ZUpper_, xU, etaQ2D);
zL = interp2(xL_old, e, ZLower_, xL, etaQ2D);

end