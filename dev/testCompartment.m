function testCompartment(Alena)
%testCompartment Tests the rigid body mass property computation for the
%compartment objects.

close all

if nargin < 0 
    warning('Need a valid instance of the Framework to run this test');
end


%% Cuboid

%Vertex coordinates
x = [ ...
    0, 4, 4, 0 ; ...
    0, 4, 4, 0 ; ...
    5, 9, 9, 5];
y = [ ...
    0, 0, 0, 0 ; ...
    5, 5, 5, 5 ; ...
    10,10,10,10];
z = [ ...
    0, 0, 3, 3 ; ...
    0, 0, 3, 3 ; ...
    0, 0, 3, 3];

%Close the surface
x = [x , x(:, 1)];
y = [y , y(:, 1)];
z = [z , z(:, 1)];

figure, grid on, hold on, box on
n = 5;
%
hAx(1) = subplot(n, 1, 1); grid on, hold on, box on
plot3(x, y, z);
title('plot3');
%
hAx(2) = subplot(n, 1, 2); grid on, hold on, box on
surf(x, y, z);
title('surf');
%
hAx(3) = subplot(n, 1, 3); grid on, hold on, box on
tri = delaunay(x(:), y(:), z(:));
TR = triangulation(tri, x(:), y(:), z(:));
[F, P] = freeBoundary(TR);
trisurf(F,P(:,1),P(:,2),P(:,3), ...
       'FaceColor','cyan','FaceAlpha',0.8);
title('deluanay -> triangulation -> freeBoundary');
%
hAx(4) = subplot(n, 1, 4); grid on, hold on, box on
TR_ = delaunayTriangulation(x(:), y(:), z(:));
[F_, P_] = freeBoundary(TR_);
trisurf(F_,P_(:,1),P_(:,2),P_(:,3), ...
       'FaceColor','blue','FaceAlpha',0.8);
title('deluanayTriangulation -> freeBoundary');
%
sf = 0.5;
hAx(5) = subplot(n, 1, 5); grid on, hold on, box on
[k, v] = boundary(x(:), y(:), z(:), sf);
trisurf(k, x(:), y(:), z(:), ...
    'FaceColor','red','FaceAlpha',0.8);
title(sprintf('boundary (shrink factor %.1f)', sf));
set(hAx, 'View', [-45, 60]);

sf = [0 : 0.1 : 1];
nSF = numel(sf);
nRow = ceil(nSF / 2);
figure
for i = 1 : nSF
    hAx(i) = subplot(nRow, 2, i); grid on, hold on, box on
    [k, v] = boundary(x(:), y(:), z(:), sf(i));
    trisurf(k, x(:), y(:), z(:), ...
        'FaceColor','red','FaceAlpha',0.8);
    title(sprintf('boundary (shrink factor %.1f)', sf(i)));
end
set(hAx, 'View', [-45, 60]);

%% Using the ALENA object

%Find the starboard wing
SW = findall(Alena, 'Name', 'Starboard Wing');

%Remove any existing compartments
deleteCompartments(SW);

%Add two compartments
addCompartment(SW, awi.model.Compartment, 'EtaLocations', [0 ; 0.5], ...
    'LocalXStart', 0.15, 'LocalXEnd', 0.65, 'PayloadFraction', 0.5);
addCompartment(SW, awi.model.Compartment, 'EtaLocations', [0.65 ; 0.95], ...
    'LocalXStart', 0.15, 'LocalXEnd', 0.65, 'PayloadFraction', 0.6);

set(SW.Compartments, 'MassDensity', 1000);
updateMassAndInertia(SW.Compartments);

%Visualise it
hF = figure;
draw(SW)
hAx = findobj(hF.Children, 'Type', 'Axes', '-not', 'Tag', '-CoordSys');
axis(hAx, 'equal');

hF(2)  = figure;
hAx(1) = axes('Parent', hF(2), 'NextPlot', 'add');
drawMesh(SW.Compartments, hAx);
axis(hAx, 'equal');

hF(3) = figure;
hAx   = axes('Parent', hF(3), 'NextPlot', 'add');
drawElement(SW.Compartments, hAx, 'DrawPayload', true);
axis(hAx, 'equal');

hF(4) = figure;
hAx   = axes('Parent', hF(4), 'NextPlot', 'add');
hg = drawMassLocations(SW.Compartments, hAx);
axis(hAx, 'equal');

end

