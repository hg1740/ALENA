function [x, z] = loadProfileDataFromFile(filename, nPoints)
%loadProfileDataFromFile Loads the aerofoil data from a text file and
%returns the (x,z) coordinates of a normalised profile.
%
% Assumptions:
%   - The x-domain is from [0 : 1] - It is assumed that the data is
%   organised by specifying
%     the upper surface (x,z) data first then the lower surface (x,z) data.
%   - The data is assumed to be in two rows seperated by a
%     single comma (',') and with one header line.
%
% TODO - Add in Airbus file format

if nargin < 2
    nPoints = 100;
end

%Validate inputs
%   - Check NACA number has been defined using a string OR a cellstr -
%   Check nPoints is a scalar integer
if ischar(filename)
    %OK
elseif iscellstr(filename)
    assert(numel(filename) == 1, 'if specified as cellstr, NACA must be single element');
    filename = filename{1};
else
    validateattributes(filename, {'string'}, {'scalartext'}, ...
        'LiftingSurface.calculateNACACoords', 'NACA');
end
validateattributes(nPoints, {'double'}, {'scalar', 'integer', ...
    'nonnan', 'finite', 'real', 'nonempty'}, ...
    'LiftingSurface.calculateNACACoords', 'nPoints');

%Open the file -> read data -> close the file
fid = fopen(filename, 'r');
assert(fid ~= -1, sprintf(['Unable to open the aerofoil ', ...
    'file ''%s''\nExiting the aerofoil import process.' ], ...
    filename));
profile_data = textscan(fid, '%f %f', 'Delimiter', ',', 'HeaderLines', 1);
fclose(fid);

%Split into (x,z)
x = profile_data{1};
z = profile_data{2};

if isempty(x) || isempty(z)
    error(i_badFormatMessage(filename));
end

%Split into upper and lower surface
index = find(or(x == 0, x == 1));
assert(numel(index) == 4, i_badFormatMessage(filename));
xU = x(index(1) : index(2));
zU = z(index(1) : index(2));
xL = x(index(3) : index(4));
zL = z(index(3) : index(4));

%Interpolate to match the desired number of points
x_c = cosspace(0, 1, nPoints);
zU  = interp1(xU, zU, x_c);
zL  = interp1(xL, zL, x_c);

%Recombine into a single vector
x = [x_c ; flipud(x_c)];
z = [zU  ; flipud(zL)];

    function msg = i_badFormatMessage(filename)
        msg = sprintf(['Aerofoil data did not meet expected '  , ...
            'format for the file ''%s''\nExiting the aerofoil ', ...
            'import process.'], filename);
    end

end