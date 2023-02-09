function [x, z] = loadCrossSectionFromFile(filename)
%loadCrossSectionFromFile Loads the cross-section (x, z) coordinates from a
%file. 
%
% Detailed Description: 
%   - The expected format of the file is two columns of data which are
%   delimited by blanks spaces or commas.
%   - Comments can be included in the file but must be denoted by #.
%   - It is expected that the first x-coordinate will be 0 or 1 and the
%     points are defined in a single continuous line around the perimeter
%     of the shape.
%
%   # I am a comment
%   # I am another comment
%   # Coord 1 | Coord 2
%       x1      y1
%       x2      y2
%       x3      y3
%       x4      y4

assert(exist(filename, 'file') == 2, sprintf(['Unable to find the file ', ...
    '''%s''. Check that the file exists and try again.']));

[x, z] = i_readDataFromFile(filename);

    function [x, z] = i_readDataFromFile(filename)
        
        fid = fopen(filename);
        strData = textscan(fid, '%s', 'Delimiter', '\n', 'CommentStyle', '#');
        strData = strData{1};
        fclose(fid);
        
        %Assume each column is delimited by blank space
        numData = cellfun(@str2num, strData, 'Unif', false);
        
        %Check format
        msg_ = sprintf(['Expected the coordinate file ''%s'' to have 2 ', ... 
            'columns of numeric data. Check file format and try again.'], filename);
        nCol = cellfun(@numel, numData);
        assert(range(nCol) == 0, msg_);
        assert(nCol(1) == 2    , msg_);        
        
        data = vertcat(numData{:});
        
        %Strip repeated points
        %[~, ind_] = unique(data, 'rows', 'stable');
        
        %if rem(numel(ind_), 2) == 0
        %    x = data(ind_, 1);
         %   z = data(ind_, 2);   
        %else
        x = data(:, 1);
        z = data(:, 2);
        %end
        
        assert(rem(numel(x), 2) == 0, ['Expected the (x, z) coordinates to ', ...
            'have an equal number of points on the upper and lower surfaces.'])
        
        %msg_ = ['Expected the x coordinates to be defined in one continuous '   ,  ...
%             'line in the range [0, 1]. At least one point should be repeated, ', ...
%             'either at x=0 or x=1.'];
        %assert(nnz(ismember(x, [0, 1]))  == 2, msg_);
        
    end

%Define coordinates in clockwise direction
if ~ispolycw(x, z)
    [x, z] = poly2cw(x, z);
end

%How is the file formatted?
msg = ['Expected the x coordinates to be defined in one continuous '   ,  ...
    'line in the range [0, 1]. At least one point should be repeated, ', ...
    'either at x=0 or x=1.'];
nEndPoint = nnz(ismember(x, [0, 1]));
%Split coordinates into upper and lower surfaces
if nEndPoint  == 2 
    if x(1) == 0
        %Starts at LE and goes clockwise
        ind1 = find(x == 1, 1);
        xU = x(1 : ind1);
        zU = z(1 : ind1);
        xL = [x(ind1 : end) ; x(1)];
        zL = [z(ind1 : end) ; z(1)];
    elseif x(1) == 1
        %Starts at TE and goes clockwise
        ind0 = find(x == 0, 1);
        xL = x(1 : ind0);
        zL = z(1 : ind0);
        xU = [x(ind0 : end) ; x(1)];
        zU = [z(ind0 : end) ; z(1)];
    else
        error('Expected the cross-section coordinates to start at x = 0 or x = 1');
    end
elseif nEndPoint == 4
    ind0 = find(x == 0);
    ind1 = find(x == 1);
    xU = x(ind0(1) : ind1(1));
    zU = z(ind0(1) : ind1(1));
    xL = x(ind0(2) : ind1(2));
    zL = z(ind0(2) : ind1(2));
else
    error(msg);
end

%Enforce cross-section to start at x = 0 and go clockwise -> Enforce unique
[xU, ind] = unique(xU);
zU        = zU(ind);
[xL, ind] = sort(xL, 'descend');
zL        = zL(ind);

%Upper and lower surfaces must be defined along the same x-positions to
%perform interpolation...
if ~isequal(xU, flipud(xL))
    zL = interp1(xL, zL, flipud(xU), 'spline');
    xL = flipud(xU);
end

x = [xU ; xL]';
z = [zU ; zL]';

end

