function [x, z] = calculateNACACoords(NACA, nPoints)
%calculateNACACoords Calculates the normalised profile coordinates of a
%generic NACA aerofoil.
%
% This function will accept the following NACA series
%   - 4-series

if nargin < 2
    nPoints = 100;
end

%Validate inputs
%   - Check NACA number has been defined using a string OR a cellstr -
%   Check nPoints is a scalar integer
if ischar(NACA)
    %OK
elseif iscellstr(NACA)
    assert(numel(NACA) == 1, 'if specified as cellstr, NACA must be single element');
    NACA = NACA{1};
else
    validateattributes(NACA, {'string'}, {'scalartext'}, ...
        'calculateNACACoords', 'NACA');
end
validateattributes(nPoints, {'double'}, {'scalar', 'integer', ...
    'nonnan', 'finite', 'real', 'nonempty'}, ...
    'calculateNACACoords', 'nPoints');

%Get NACA number only - Need it as a char array
nacaNum = char(strrep(NACA, 'NACA',''));

%             %Define function handle f =
%             str2func(sprintf('i_naca%iseries', numel(nacaNum)));
%             %Calculate the profile profile = f(nacaNum, nPoints);

%Determine NACA series
series = i_whatSeries(nacaNum);

%Define the vector of x-points using cosine spacing
x_c = cosspace(0, 1, nPoints);

%Define function handle
switch series
    case '4-series'
        [x, z] = i_naca4series(nacaNum, x_c);
    otherwise
        error('Uknown NACA number encountered');
end

    function series = i_whatSeries(nacaNum)
        %i_whatSeries Determines the NACA aerofoil series based on the NACA
        %number.
        %
        % numel  --    series
        %   4        '4-series' 5        '5-series'
        
        if numel(nacaNum) == 4
            series = '4-series';
        else
            error('Unknown NACA number encountered');
        end
    end

    function [x, z] = i_naca4series(nacaNum, x_c)
        %i_naca4series Calculates the (x,z) profile coordinates for a 4
        %digit NACA aerofoil at the non-dimensionsal chordwise position
        %'x_c'.
        
        %Derive NACA terms
        m = str2double(nacaNum(1))   ./ 100;
        p = str2double(nacaNum(2))   ./ 10;
        t = str2double(nacaNum(3:4)) ./ 100;
        
        %Calculate the thickness line 
        a0 = 0.2969;
        a1 = -0.1260;
        a2 = -0.3516;
        a3 = 0.2843;
        a4 = -0.1036;
        zt = 5 .* t .* (a0 .* (x_c).^(0.5) + a1 .* x_c + a2 .* x_c.^2 + a3 .* x_c.^3 + a4 .* x_c.^4);
        zt(end) = 0;    % ensure zero trailing edge!
        
        %Logical indexing
        idx_1 = (x_c <= p); % 0 <= x <= p
        idx_2 = (x_c >= p); % p <= x <= 1
                
        %Calculate the y-coordinates of the mean camber line
        zc(idx_1, 1) = (m/p^2) .* (2 .* p .* x_c(idx_1) - x_c(idx_1).^2);
        zc(idx_2, 1) = (m./((1-p).^2)) .* ((1-2.*p) + 2.*p .* x_c(idx_2) - x_c(idx_2).^2);
        
        %Calculate the rate of change of the camber line w.r.t x
        dzc_dx(idx_1, 1) = (2 .* m ./ p.^2).*(p-x_c(idx_1));
        dzc_dx(idx_2, 1) = (2 .* m ./ ((1 - p).^2)).*(p-x_c(idx_2));
        
        %Calculate the angle of the mean camber line
        theta = atan(dzc_dx);
        
        %Define upper and lower surfaces
        xU = x_c - zt .* sin(theta);    %Upper-x
        xL = x_c + zt .* sin(theta);    %Lower-x
        zU = zc + zt .* cos(theta);     %Upper-z
        zL = zc - zt .* cos(theta);     %Lower-z                
        
        %Combine upper and lower surfaces into a single vector
        x = [xU ; flipud(xL)];
        z = [zU ; flipud(zL)];
        
        %Sometimes the 'xU' and 'xL' can be less than zero due to sin/cos
        %terms - set to abs
        x = abs(x);
        
    end
end