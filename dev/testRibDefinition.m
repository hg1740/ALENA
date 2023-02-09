function testRibDefinition(LS, etaRib, lambdaRib)
%testRibDefition Tests the definition of Rib objects using an arbitrary
%volume slice.

close all
   
debugMode = true; %toggle for showing helpful plots
    
if nargin < 2
    etaRib = [0.01, 0.1, 0.4, 0.6, 0.8];
end
if nargin < 3
    lambdaRib = 40 * rand(size(etaRib));
end

validateattributes(LS, {'awi.model.LiftingSurface'}, {'scalar'}, 'testRibDefinition', 'LS');

Spars = find(LS.Children, 'Type', 'Spar');
if isempty(Spars)
    warning('Unable to calculate Ribs without Spar objects.');
    return
end

%Strip any ribs defined at eta = 0 and eta = 1, these are added manually
idx       = ~ismember(etaRib, [0, 1]);
etaRib    = etaRib(idx);
lambdaRib = lambdaRib(idx);

%Get the structural layout
[Ribs, Stringers] = defineStructuralLayout(LS, Spars);

if debugMode
    hF  = figure('Name', 'Spar geometry and rib intersection planes');
    hAx = axes('Parent', hF, 'NextPlot', 'add');
    drawElement(Spars, hAx);
end

%Check all 'Rib' objects have the same number of points
assert(all(range(arrayfun(@(R) numel(R.X), Ribs)) == 0), ['Expected the ', ...
    'rib CrossSection objects to have the same number of points.']);

%Define the intersection planes
nRib    = numel(etaRib);
[p0, n] = defineRibIntersectionPlanes(LS, Ribs, etaRib, lambdaRib, nRib);

%% Get coordinates of control Cross-Sections in global frame

[x, y, z] = calculateGlobalCoords(Ribs);

nControlSection = size(Ribs, 2);

%Combine rib coordinates at each spanwise station
sz = size(Ribs);
x = reshape(num2cell(x, 2), sz);
y = reshape(num2cell(y, 2), sz);
z = reshape(num2cell(z, 2), sz);
%   - Split into upper and lower surfaces (indexing is done to avoid
%   duplicate data points)
nP_ = (numel(x{1}) - 1) / 2;
xU  = cellfun(@(v) v(2 : nP_), x, 'Unif', false);
yU  = cellfun(@(v) v(2 : nP_), y, 'Unif', false);
zU  = cellfun(@(v) v(2 : nP_), z, 'Unif', false);
xL  = flipud(cellfun(@(v) v(nP_ + 2 : end - 1), x, 'Unif', false));
yL  = flipud(cellfun(@(v) v(nP_ + 2 : end - 1), y, 'Unif', false));
zL  = flipud(cellfun(@(v) v(nP_ + 2 : end - 1), z, 'Unif', false));
%   - Combine
x = arrayfun(@(i) horzcat(x{1, i}(1), xU{:, i}, x{1, i}(nP_ + 1), xL{:, i}, x{1, i}(1)), 1 : nControlSection, 'Unif', false);
y = arrayfun(@(i) horzcat(y{1, i}(1), yU{:, i}, y{1, i}(nP_ + 1), yL{:, i}, y{1, i}(1)), 1 : nControlSection, 'Unif', false);
z = arrayfun(@(i) horzcat(z{1, i}(1), zU{:, i}, z{1, i}(nP_ + 1), zL{:, i}, z{1, i}(1)), 1 : nControlSection, 'Unif', false);
x = vertcat(x{:});
y = vertcat(y{:});
z = vertcat(z{:});

nP = size(x, 2);

%% Calculate rib coordinates using vector-plane intersection

%Calculate direction vector of each cross-section vertex along the span
[r0, r_mod, r_norm] = getRibVector(x, y, z);

%Calculate intersection coordinates for each rib
[sI, xI, yI, zI] = vectorPlaneIntersection(n, p0, r_norm, r0);

%Intersection occurs between [0 : 1] of the vector length
dsI = sI ./ r_mod;
idx = and(dsI <= 1, dsI >= 0);

%Stash coordinates
numCoords = nnz(idx) / nP;
xRib = reshape(xI(idx), [nP, numCoords]);
yRib = reshape(yI(idx), [nP, numCoords]);
zRib = reshape(zI(idx), [nP, numCoords]);

%Add the coordinates at eta = [0, 1]
xRib = [x(1, :)' , xRib , x(end, :)'];
yRib = [y(1, :)' , yRib , y(end, :)'];
zRib = [z(1, :)' , zRib , z(end, :)'];

if debugMode
    plot3(hAx, xRib, yRib, zRib);
end

%% Update Spar coordinates to include intersection points with ribs

%Get spar coordinates (upper and lower surface)
[xSp, ySp, zSp] = calculateSparCoords(Spars);

if debugMode
    hF  = figure('Name', 'Spar control points and rib intersection points');
    hAx = axes('Parent', hF, 'NextPlot', 'add');
    hold on
    plot3(hAx, xRib, yRib, zRib);
    plot3(hAx, xSp{1}', ySp{1}', zSp{1}', 'r-', 'Marker', 'o', 'MarkerFaceColor', 'r', 'MarkerEdgeColor', 'k');
    plot3(hAx, xSp{2}', ySp{2}', zSp{2}', 'b-', 'Marker', 'o', 'MarkerFaceColor', 'b', 'MarkerEdgeColor', 'k');
end

xSp = vertcat(xSp{:})';
ySp = vertcat(ySp{:})';
zSp = vertcat(zSp{:})';

%Calculate direction vector of each cross-section vertex along the span
[r0, r_mod, r_norm] = getRibVector(xSp, ySp, zSp);

%Calculate intersection coordinates for each rib
[sI, xI, yI, zI] = vectorPlaneIntersection(n, p0, r_norm, r0);

%Intersection occurs between [0 : 1] of the vector length
dsI = sI ./ r_mod;
idx = and(dsI <= 1, dsI >= 0);

if debugMode
    plot3(hAx, xI(idx), yI(idx), zI(idx), 'Marker', '^', 'MarkerFaceColor', 'g', 'LineStyle', 'none', 'MarkerEdgeColor', 'k')
end

nSp   = numel(Spars);
nSpar = 2 * nSp;
xRibSpar = reshape(xI(idx), [nSpar, nRib]);
yRibSpar = reshape(yI(idx), [nSpar, nRib]);
zRibSpar = reshape(zI(idx), [nSpar, nRib]);

%Add intersection points and sort in spanwise direction
%   - TODO: Add switch statement for 'LS.SpanVector'
xSp = [xSp ; xRibSpar'];
ySp = [ySp ; yRibSpar'];
zSp = [zSp ; zRibSpar'];
[~, ind] = sort(ySp(:, 1), 'ascend');
xSp = xSp(ind, :)';
ySp = ySp(ind, :)';
zSp = zSp(ind, :)';

%Append spar coordinates at eta = [0, 1] to the rib-spar intersection set
xRibSpar = [xSp(:, 1), xRibSpar, xSp(:, end)];
yRibSpar = [ySp(:, 1), yRibSpar, ySp(:, end)];
zRibSpar = [zSp(:, 1), zRibSpar, zSp(:, end)];

%Set up coordinates in format required for mesh seed generation
xSp = cat(3, xSp(1 : nSp, :), xSp(nSp + 1 : end, :));
ySp = cat(3, ySp(1 : nSp, :), ySp(nSp + 1 : end, :));
zSp = cat(3, zSp(1 : nSp, :), zSp(nSp + 1 : end, :));

%% Add intersection points to the rib coordinates

%Split into upper and lower portions
xUrs = xRibSpar(1 : nSp, :);
yUrs = yRibSpar(1 : nSp, :);
zUrs = zRibSpar(1 : nSp, :);
xLrs = xRibSpar(nSp + 1 : end, :);
yLrs = yRibSpar(nSp + 1 : end, :);
zLrs = zRibSpar(nSp + 1 : end, :);

%Append the rib-spar intersection point - first and last spar already
%accounted for
nPoints = (nP - 1) / 2;
xRibU   = [xRib(1 : nPoints, :) ; xUrs(2 : end - 1, :)];
yRibU   = [yRib(1 : nPoints, :) ; yUrs(2 : end - 1, :)];
zRibU   = [zRib(1 : nPoints, :) ; zUrs(2 : end - 1, :)];
% xRibL   = [xRib(nPoints + 1 : end - 1, :) ; xLrs(2 : end - 1, :)];
% yRibL   = [yRib(nPoints + 1 : end - 1, :) ; yLrs(2 : end - 1, :)];
% zRibL   = [zRib(nPoints + 1 : end - 1, :) ; zLrs(2 : end - 1, :)];

xRibL   = [flipud(xRib(nPoints + 1 : end - 1, :)) ; xLrs(2 : end - 1, :)];
yRibL   = [flipud(yRib(nPoints + 1 : end - 1, :)) ; yLrs(2 : end - 1, :)];
zRibL   = [flipud(zRib(nPoints + 1 : end - 1, :)) ; zLrs(2 : end - 1, :)];

%Check for duplicates and sort in the streamwise (global-X) direction
%[~, ia, ib] = unique([xRibU(:), yRibU(:), zRibU(:)], 'rows');
%[xRibU, ia, ~] = uniquetol(xRibU, 1e-4, 'ByRows', true);
[xRibU, yRibU, zRibU] = i_uniqueRibPoints(xRibU, yRibU, zRibU);
[xRibL, yRibL, zRibL] = i_uniqueRibPoints(xRibL, yRibL, zRibL);

    function [ux, uy, uz] = i_uniqueRibPoints(x, y, z)
         
        [~, iax, ~] = unique(x, 'rows');
        %[~, iax, ~] = unique(x, 'rows', 'stable');
        %[~, iay, ~] = unique(y, 'rows', 'stable');
        %[~, iaz, ~] = unique(z, 'rows', 'stable');
        
        %assert(~any(range([iax, iay, iaz], 2)), 'Unexpected sort in the GFEM meshing code.');
        
        ux = x(iax, :);
        uy = y(iax, :);
        uz = z(iax, :);
        
        %[ux, ind_] = sort(ux(:, 1), 'ascend');
        %uy = uy(ind_, :);
        %uz = uz(ind_, :);
        
    end

%[xRibL, ia, ~] = uniquetol(xRibL, 1e-4, 'ByRows', true);
%[xRibL, ia, ~] = unique(xRibL, 'rows');
%xRibL = flipud(xRibL);
%yRibL = flipud(yRibL);
%zRibL = flipud(zRibL);

%% Define spar mesh

nEta = size(xSp, 2);
%Mesh size
sw = LS.ShellWidth;
ar = LS.ShellAR;
if isempty(sw)
    sw = 0.1;
end
if isempty(ar)
    ar = 1;
end
sh = ar * sw;
type = 'uniform';
[meshCoords, nShellHeight] = i_getVerticalMeshSeed(xSp, ySp, zSp, nSp, nEta, sh, type);

sparMesh = i_getSparMesh(xSp, ySp, zSp, meshCoords, nSp, nEta, nShellHeight, sw, type);

%Plot to check?
if debugMode
    coords = vertcat(sparMesh{:});
    i_flatten = @(x) reshape(x, [numel(x), 1]);
    x = i_flatten(coords(:, :, 1));
    y = i_flatten(coords(:, :, 2));
    z = i_flatten(coords(:, :, 3));
    hF = figure('Name', 'Spar Spanwise Mesh');
    hAx = axes('Parent', hF, 'NextPlot', 'add');
    xlabel(hAx, 'X [m]');
    ylabel(hAx, 'Y [m]');
    zlabel(hAx, 'Z [m]');
    %hg = drawElement(Spars, hAx);
    %set(hg, 'FaceAlpha', 0.1);
    hg(1) = plot3(hAx, x, y, z, ...
        'Marker'         , '^', ...
        'MarkerFaceColor', 'g', ...
        'MarkerEdgeColor', 'k', ...
        'LineStyle'      , 'none', 'DisplayName', 'Spar Nodes');
    legend(hAx, hg, get(hg, {'DisplayName'}));
end

%% Define rib mesh

%Find index number of rib-spar intersections
indU = cell(3, nSp - 2);
indL = indU;
for iSp = 1 : nSp - 2
    indU{1, iSp} = find(all((xRibU - xUrs(iSp + 1, :)) == 0, 2));
    indU{2, iSp} = find(all((yRibU - yUrs(iSp + 1, :)) == 0, 2));
    indU{3, iSp} = find(all((zRibU - zUrs(iSp + 1, :)) == 0, 2));
    indL{1, iSp} = find(all((xRibL - xLrs(iSp + 1, :)) == 0, 2));
    indL{2, iSp} = find(all((yRibL - yLrs(iSp + 1, :)) == 0, 2));
    indL{3, iSp} = find(all((zRibL - zLrs(iSp + 1, :)) == 0, 2));    
end
assert(all(range(cellfun(@numel, indU)) == 0), ['Expected all ribs to have ', ...
    'the same number of spar intersections.']);
assert(all(range(cellfun(@numel, indL)) == 0), ['Expected all ribs to have ', ...
    'the same number of spar intersections.']);
indU = cell2mat(indU);
indL = cell2mat(indL);
assert(all(arrayfun(@(i) range(indU(:, i)), 1 : nSp - 2) == 0), ['Expected the rib-spar intersections to ', ...
    'have the same index for (x,y,z).']);
assert(all(arrayfun(@(i) range(indU(:, i)), 1 : nSp - 2) == 0), ['Expected the rib-spar intersections to ', ...
    'have the same index for (x,y,z).']);
indU = indU(1, :);
indL = indL(1, :);
lb = [1, indU];
ub = [indU, size(xRibU, 1)];

%Determine node positions on rib perimeter
for iBay = 1 : nSp - 1
    
    %Split rib perimeter at spar intersection points
    xU_r = xRibU(lb(iBay) : ub(iBay), :);
    yU_r = yRibU(lb(iBay) : ub(iBay), :);
    zU_r = zRibU(lb(iBay) : ub(iBay), :);
    xL_r = xRibL(lb(iBay) : ub(iBay), :);
    yL_r = yRibL(lb(iBay) : ub(iBay), :);
    zL_r = zRibL(lb(iBay) : ub(iBay), :);
    
    %Calculate straight-line distance along the perimeter
    [~, r_modU, ~] = getRibVector(xU_r, yU_r, zU_r);
    [~, r_modl, ~] = getRibVector(xL_r, yL_r, zL_r);
    rU = [zeros(1, nRib +  2) ; cumsum(r_modU)];
    rL = [zeros(1, nRib +  2) ; cumsum(r_modl)];
    
    %Determine shell width from length of rib bay perimeter
    rMax     = max([rU(end, :), rL(end, :)]);
    nRibNode = ceil(rMax / sw);
    
    %Interpolate the rib perimeter at the node positions    
    coordsU = zeros(nRibNode - 2, size(rU, 2), 3);
    coordsL = zeros(nRibNode - 2, size(rU, 2), 3);
    %xU_r_n = zeros(nRibNode - 2, size(rU, 2));
    %yU_r_n = xU_r_n;
    %zU_r_n = xU_r_n;
    %xL_r_n = xU_r_n;
    %yL_r_n = xU_r_n;
    %zL_r_n = xU_r_n;
    for iRib = 1 : size(rU, 2)
        r_interp_U = linspace(0, 1, nRibNode) * rU(end, iRib);
        r_interp_L = linspace(0, 1, nRibNode) * rL(end, iRib);
        ribNodeCoord_U = interp1(rU(:, iRib), [xU_r(:, iRib), yU_r(:, iRib), zU_r(:, iRib)], r_interp_U(2 : end - 1)', 'linear', 'extrap'); 
        ribNodeCoord_L = interp1(rL(:, iRib), [xL_r(:, iRib), yL_r(:, iRib), zL_r(:, iRib)], r_interp_L(2 : end - 1)', 'linear', 'extrap'); 
%         xU_r_n(:, iRib) = ribNodeCoord_U(:, 1);
%         yU_r_n(:, iRib) = ribNodeCoord_U(:, 2);
%         zU_r_n(:, iRib) = ribNodeCoord_U(:, 3);
%         xL_r_n(:, iRib) = ribNodeCoord_L(:, 1);
%         yL_r_n(:, iRib) = ribNodeCoord_L(:, 2);
%         zL_r_n(:, iRib) = ribNodeCoord_L(:, 3);
        coordsU(:, iRib, :) = permute(ribNodeCoord_U, [1, 3, 2]); 
        coordsL(:, iRib, :) = permute(ribNodeCoord_L, [1, 3, 2]); 
    end
    
    %Calculate node positions for rib face mesh
    
end


%Resample the rib perimeter at the desired node positions

end

%% Vector-plane intersection methods

function [p0, n] = defineRibIntersectionPlanes(LS, Ribs, etaRib, lambdaRib, nRib)
%defineRibIntersectionPlanes Defines the rib intersection planes in terms
%of a normal vector and a point for each rib plane.

%Get the bounding box of the rib coordinate data
%   - Make the bounding box slightly bigger to account for rotations of the
%     intersection plane.
[x, y, z] = calculateGlobalCoords(Ribs);
ribCoords = horzcat(x(:), y(:), z(:));
minData   = 1.25 .* min(ribCoords, [], 1);
maxData   = 1.25 .* max(ribCoords, [], 1);

%Define intersection planes - TODO: Add checkfor SpanVector
%   - Use corners of bounding box
xd = repmat([minData(1) ; maxData(1) ; maxData(1) ; minData(1)], [1, nRib]);
yd = repmat(etaRib * LS.Span, [4, 1]);
zd = repmat([minData(3) ; minData(3) ; maxData(3) ; maxData(3)], [1, nRib]);
%   - Get centre of coordinates (so we can
xc = mean(xd, 1);
yc = mean(yd, 1);
zc = mean(zd, 1);
%   - Rotate by 'lambdaRib' about the centre of each plane
for iRib = 1 : nRib
    rot = [ ...
        cosd(lambdaRib(iRib)) , sind(lambdaRib(iRib)), 0  ; ...
        -sind(lambdaRib(iRib)), cosd(lambdaRib(iRib)), 0  ; ...
        0          ,           0          , 1 ];
    coords = [xd(:, iRib), yd(:, iRib), zd(:, iRib)]' - [xc(iRib) ; yc(iRib) ; zc(iRib)];
    transf_coords = rot * coords;
    xd(:, iRib) = transf_coords(1, :) + xc(iRib);
    yd(:, iRib) = transf_coords(2, :) + yc(iRib);
    zd(:, iRib) = transf_coords(3, :) + zc(iRib);
end
%   - Plot it
fcs = [ ...
    1 : 4 : 4 * nRib ; ...
    2 : 4 : 4 * nRib ; ...
    3 : 4 : 4 * nRib ; ...
    4 : 4 : 4 * nRib]';
patch(gca, 'Vertices', [xd(:), yd(:), zd(:)], 'Faces', fcs, 'FaceColor', 'b', 'FaceAlpha', 0.1);

%Define planes in terms of a normal vector and a point in the plane
r1 = [xd(2, :) - xd(1, :) ; yd(2, :) - yd(1, :) ; zd(2, :) - zd(1, :)];
r2 = [xd(4, :) - xd(1, :) ; yd(4, :) - yd(1, :) ; zd(4, :) - zd(1, :)];
n  = cross(r2, r1);
n  = n ./ vecnorm(n);
p0 = [xd(1, :) ; yd(1, :) ; zd(1, :)];

end

function [r0, r_mod, r_norm] = getRibVector(x, y, z)

%Construct the vector that runs along the lofted cross-sections
dX     = diff(x, 1);
dY     = diff(y, 1);
dZ     = diff(z, 1);
r_mod  = sqrt(sum(cat(3, dX.^2, dY.^2, dZ.^2), 3));
r_norm = cat(3, dX ./ r_mod, dY ./ r_mod, dZ ./ r_mod);
r0     = cat(3, x(1 : end - 1, :), y(1 : end - 1, :), z(1 : end - 1, :));

end

function [sI, xI, yI, zI] = vectorPlaneIntersection(n, p0, u, r0)

%Split into (x,y,z) components to allow vectorised operations
nx  = permute(n(1, :) , [1, 3, 2]);
ny  = permute(n(2, :) , [1, 3, 2]);
nz  = permute(n(3, :) , [1, 3, 2]);
p0x = permute(p0(1, :), [1, 3, 2]);
p0y = permute(p0(2, :), [1, 3, 2]);
p0z = permute(p0(3, :), [1, 3, 2]);
r0x = r0(:, :, 1);
r0y = r0(:, :, 2);
r0z = r0(:, :, 3);
ux  = u(:, :, 1);
uy  = u(:, :, 2);
uz  = u(:, :, 3);

%W = r0 - P0. We want dot(n, -w) so we set w = -(r0 - P0) = (P0 - r0)
wx = p0x - r0x;
wy = p0y - r0y;
wz = p0z - r0z;

%Calculate fraction of vector 'u' which the intersection occurs at (sI)
n_dot_w = wx .* nx + wy .* ny + wz .* nz;
n_dot_r = ux .* nx + uy .* ny + uz .* nz;
sI = n_dot_w ./ n_dot_r;

%Calculate intersection coordinates
xI = r0x + sI .* ux;
yI = r0y + sI .* uy;
zI = r0z + sI .* uz;

end

%% Defining mesh seed

function [meshCoords, nShellHeight] = i_getVerticalMeshSeed(xSp, ySp, zSp, nSp, nEta, sh, type)
%i_getVerticalMeshSeed Returns the (x,y,z) coordinates of
%the grid nodes at each spanwise 'eta' position for each
%spar. The position of the grid nodes is calculated based
%on the desired shell dimensions or the number of shells.

%Determine the number shells required based on the spar height
%at each 'eta' position
%   - 'r' is the vector between the upper and lower surfaces of
%     the spar at every vertex along the span
%   - 'r_mod' is the straight line distance between the upper
%     and lower vertices at every location along the span.
%   - 'r_norm' is the direction vector
[~, r_mod, r_norm] = i_getVector(xSp, ySp, zSp, 3, 3);
nShellHeight = max(ceil(r_mod ./ sh), [], 2);

if strcmp(type, 'uniform')
    nShellHeight = repmat(max(nShellHeight), size(nShellHeight));
end

%Mesh each spar individually...

%Mesh seed along the height of the spar
%   - meshCoords = [nShell, nEta, 3] Dim 3 => (x, y, z)
meshCoords = arrayfun(@(n) zeros(n + 1, nEta, 3), nShellHeight, 'Unif', false);
for iSp = 1 : nSp
    %Define seed for new points - start at upper vertices
    r0 = [xSp(iSp, :, 1) ; ySp(iSp, :, 1) ; zSp(iSp, :, 1)];
    ds = cumsum(repmat(r_mod(iSp, :) ./ nShellHeight(iSp), [nShellHeight(iSp) - 1, 1]), 1);
    %Mesh uses lower/upper vertices as initial seed
    meshCoords{iSp}(1  , :, 1) = xSp(iSp, :, 1);
    meshCoords{iSp}(1  , :, 2) = ySp(iSp, :, 1);
    meshCoords{iSp}(1  , :, 3) = zSp(iSp, :, 1);
    meshCoords{iSp}(end, :, 1) = xSp(iSp, :, 2);
    meshCoords{iSp}(end, :, 2) = ySp(iSp, :, 2);
    meshCoords{iSp}(end, :, 3) = zSp(iSp, :, 2);
    %Intermediate coords found by using spar direction vector
    %   - r_i = r_0 + r_norm * ds
    meshCoords{iSp}(2 : end - 1, :, 1) = r0(1, :) + r_norm(iSp, :, 1) .* ds;
    meshCoords{iSp}(2 : end - 1, :, 2) = r0(2, :) + r_norm(iSp, :, 2) .* ds;
    meshCoords{iSp}(2 : end - 1, :, 3) = r0(3, :) + r_norm(iSp, :, 3) .* ds;
end

end

function mesh = i_getSparMesh(xSp, ySp, zSp, meshSeed, nSp, nEta, nShellHeight, sw, type)
%i_getSparMesh Returns the (x,y,z) coordinats of the spar
%mesh based on the mesh seed and the shell dimensions

%How many shells are required along the span?
[~, r_mod, ~] = i_getVector(xSp, ySp, zSp, 2, 4);
nShell        = max(ceil(r_mod ./ sw), [], 3);
if strcmp(type, 'uniform')
    nShell = repmat(max(nShell, [], 1), [size(nShell, 1), 1]);
end
nShellSpan = sum(nShell, 2);

%Calculate mesh points on a bay-by-bay basis for each spar
mesh = arrayfun(@(i) zeros(nShellHeight(i) + 1, nShellSpan(i) + 1, 3), 1 : nSp, 'Unif', false);
for iSp = 1 : nSp
    %Direction vector between each spanwise point
    r      = diff(meshSeed{iSp}, [], 2);
    r_mod  = vecnorm(r, 2, 3);
    r_norm = r ./ r_mod;
    %How many elements per bay?
    %                     nShellPerBay = max(ceil(r_mod ./ sw), [], 1);
    nShellPerBay = nShell(iSp, :);
    ub      = cumsum(nShellPerBay);
    ub(end) = ub(end);
    lb      = [1, ub(1 : end - 1) + 1];
    %Loop through bays and generate new points
    for iE = 1 : nEta - 1
        r0   = squeeze(meshSeed{iSp}(:, iE, :));
        r_n_ = squeeze(r_norm(:, iE, :));
        ds   = cumsum([zeros(nShellHeight(iSp) + 1, 1), repmat(r_mod(:, iE) ./ nShellPerBay(:, iE), [1, nShellPerBay(iE) - 1])], 2);
        mesh{iSp}(:, lb(iE) : ub(iE), 1) = r0(:, 1) + r_n_(:, 1) .* ds;
        mesh{iSp}(:, lb(iE) : ub(iE), 2) = r0(:, 2) + r_n_(:, 2) .* ds;
        mesh{iSp}(:, lb(iE) : ub(iE), 3) = r0(:, 3) + r_n_(:, 3) .* ds;
    end
    mesh{iSp}(:, end, :) = meshSeed{iSp}(:, end, :);
end

end

function [r, r_mod, r_norm] = i_getVector(x, y, z, dimDiff, dimCat)
%i_getVector Calculates the vector ('r') between a set of
%points. Also returns the length of the straight line
%('r_mod') and the normalised direction vector ('r_norm').
r      = cat(dimCat, diff(x, [], dimDiff), diff(y, [], dimDiff), diff(z, [], dimDiff));
r_mod  = vecnorm(r, 2, dimCat);
r_norm = r ./ r_mod;
end

