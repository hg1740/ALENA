function crossSectionGeometryGenerator(coords)
%crossSectionGeometryGenerator Uses NURBS surfaces to generate the geometry
%of a cross-section.
%
%
% References:
%   [1]. "A software tool for generic parameterized aircraft design", 2006
%        Sotirios S. Sarakinos 1, Ioannis M. Valakos 1, Ioannis K. Nikolos.
%   [2]. "The NURBS book", 2nd Edition, 

%Check we are running 2016b or higher
matlabVersion = version('-release');
releaseNum    =  sscanf(matlabVersion, '%f');
if releaseNum < 2016
    error('Function requires MATLAB release 2016b or higher');
elseif releaseNum == 2016 && strcmp(matlabVersion, '2016a')
    error('Function requires MATLAB release 2016b or higher');
end

close all

%Set up graphics objects
hF = figure('Name', 'Define control points for Bezier curve', ...
    'Visible'              , 'off'           , ...
    'WindowButtonDownFcn'  , @cbMouseClicked , ...
    'WindowButtonUpFcn'    , @cbMouseReleased, ...
    'WindowButtonMotionFcn', @cbMouseMotion);
hAx = axes('Parent', hF, ...
    'NextPlot', 'add'  , ...
    'XLim'    , [0 100], ...
    'YLim'    , [0 100]);
setappdata(hF, 'MouseDepressedOnMarker', false);

%% Get points for Bezier curve
if nargin < 1
    answer = questdlg('Generate points through GUI or use defaults?', ...
        'Bezier Curve Setup', ...
        'Use the GUI', 'Use defaults', 'Use the GUI');
    if isempty(answer) %Escape route
        return
    end
    switch answer
        case 'Use the GUI'
            %Hijack the built-in 'ginput' function and extract points from
            %user input
            xlabel(hAx, 'X [m]'); 
            ylabel(hAx, 'Y [m]'); 
            zlabel(hAx, 'Z [m]');
            [cpX, cpY, cpZ] = custom_ginput;
            coords = [cpX, cpY, cpZ]';
        case 'Use defaults'
            %Use 2D coordinates as default
            coords = [ ...
                10, 30, 60, 90 ; ...
                10, 40, 30, 15 ; ...
                0, 0, 0, 0];
    end
end
hF.Visible = 'on';

%Parse
validateattributes(coords, {'numeric'}, {'2d', 'nrows', 3, 'nonempty', ...
    'nonnan', 'finite', 'real'}, 'crossSectionGeometryGenerator', 'coords');

%Stash the Bezier curve control points in the figure appdata
setappdata(hF, 'BezierControlPointCoords', coords');

%% Bezier curves

%How many control points?
nPower = size(coords, 2);

%Parametric coordinates
nPoints = 50;
u = linspace(0, 1, nPoints);

%Calculate Bernstein polynomials
Bi = getBernsteinPoly(nPower, u);

%Stash the Bernstein polynomials
setappdata(hF, 'BernsteinPolynomials', Bi);

%Recover the Bezier curve coordinates
bezierCurve = coords * Bi;

%% Plot the data
hCP = findobj(hAx, 'Type', 'line', '-and', 'Tag', 'BezierControlPoint');
if isempty(hCP)
    hCP  = plot3(hAx, coords(1, :), coords(2, :), coords(3, :), ...
        'Marker'         , 'o', ...
        'MarkerFaceColor', 'b', ...
        'MarkerEdgeColor', 'k', ...
        'LineStyle'      , '-'); 
end
hBez = plot3(hAx, bezierCurve(1, :), bezierCurve(2, :), bezierCurve(3, :), ...
    'LineStyle', '-', ...
    'Color'    , [1, 0 0]);
hCS = patch(hAx, ...
    [bezierCurve(1, :), bezierCurve(1, 1)], ...
    [bezierCurve(2, :), bezierCurve(2, 1)], ...
    [bezierCurve(3, :), bezierCurve(3, 1)], ...
    'FaceColor', 'red', ...
    'FaceAlpha', 0.1  , ...
    'LineStyle', 'none');
%Update figure name and add helpful annotations
set(hF, 'Name', 'Customisable Bezier Curve');
title(hAx, 'Click and drag the control points to change the Bezier curve');
legend(hAx, [hCP, hBez, hCS], {'Control Points', 'Bezier Curve', 'Watertight Cross-section'});

%Stash the handle of the Bezier curve control points in the figure appdata
setappdata(hF, 'BezierControlPointHandle', hCP);
setappdata(hF, 'BezierCurveHandle'       , hBez);
setappdata(hF, 'CrossSectionPatchHandle' , hCS);

return

%Plot
hF  = figure('Name', 'Bezier Curves');
%Bernstein Polynomials
hAx(1) = axes('Parent', hF, 'NextPlot', 'add', 'OuterPosition', [0  , 0, 0.5, 1]);
hL  = plot(u, Bi);
str = arrayfun(@(n) sprintf('n = %i', n), pow, 'Unif', false);
legend(hAx(1), hL, str, 'Location', 'NorthEast');
xlabel(hAx(1), 'u [-]');
ylabel(hAx(1), 'B_{i} [-]');
title(hAx(1), sprintf('Bernestien polynomials in increasing power up to n = %i', nPower));
%Bezier curves
hAx(2) = axes('Parent', hF, 'NextPlot', 'add', 'OuterPosition', [0.5, 0, 0.5, 1]);
hL(1) = plot(hAx(2), coords(1, :), coords(2, :), 'Marker', 'o', 'MarkerFaceColor', 'b', 'MarkerEdgeColor', 'k', 'LineStyle', 'none');
hL(2) = plot(hAx(2), bezierCurve(1, :), bezierCurve(2, :), 'b-');
xlabel(hAx(2), 'X [m]');
ylabel(hAx(2), 'Y [-]');
title(hAx(2), sprintf('Bezier curve for n = %i control points', nPower));
axis(hAx(2), 'equal');

%% Inputs

% return

%Parameters
nSections = 2;
Lref = 10;

%Non-dimensional parameters for controlling the reference volume (rv) size
lrv_ = 0.5; %Length
hrv_ = 0.1; %Height
trv_ = 0.3; %Width

%(x,y,z) locations of local coordinate systems for cross-sections
xrv_ = [0, 1];
yrv_ = [0, 1];
zrv_ = [0, 1];

%% Dependent parameters

%Physical dimensions of the reference volume
lrv = Lref * lrv_;  %[m], length
hrv = hrv_ * lrv;   %[m], half height
trv = trv_ * hrv;   %[m], half width

%% Plotting

%Set up graphics objects
hF  = figure('Name', 'Reference Volume');
hAx = axes('Parent', hF, 'NextPlot', 'add');
xlabel(hAx, 'X [m]');
ylabel(hAx, 'Y [m]');
zlabel(hAx, 'Z [m]');
axis(hAx, 'equal');
view(hAx, [-40, 50]);
hAx.Box = 'on';

%Plotting
drawRefVolume(lrv, hrv, trv, hAx);

end

%Bernstein Polynomials
function Bi = getBernsteinPoly(n, u)
%getBernsteinPoly Calculates the Bernstein polynomials of power 'n' using
%parametric coordinates 'u'.

%Calculate Bernstein Polynomials
% B_(i,n) = (n! / (i! * (n - i!))) * u^i * (1 - u^i) ^ (n - i)
%   - Binomial coefficients
pow_ = (0 : n - 1)';
Ki = factorial(n - 1) ./ (factorial(pow_) .* factorial(n - 1 - pow_));
%   - 'u' raised to powers of i
pow    = (1 : n)';
uI  = arrayfun(@(i) u .^ i, pow - 1, 'Unif', false);
uI  = vertcat(uI{:});
% uI = u .^ (pow - 1);
%   - (1 - u) raise to powers of (n - i)
uNI = arrayfun(@(ni) (1 - u) .^ ni, n - pow, 'Unif', false);
uNI = vertcat(uNI{:});
% uNI = (1 - u) .^ (nPower - pow);
%   - Combine to get the actual Bernstein polynomials
Bi = Ki .* uI .* uNI;

end

%Plotting
function hg = drawRefVolume(lrv, hrv, trv, hAx, origin)
%drawRefVolume Draws the reference volume parallelipiped.
%
% Inputs:
%   * 'lrv'    - Length of the parallelipiped.
%   * 'hrv'    - Half height of the parallelipiped.
%   * 'trv'    - Half width of the parallelipiped.
%   * 'hAx'    - Handle to MATLAB axes object for plotting.
%   * 'origin' - Origin of the parallelipiped.

if nargin < 5
    origin = [0 0 0];
end

clr  = 'red';   % face color
alfa = 0.1;     % face transparency

%Generate vertex coordinate data
%   - Points are ordered starting in the positive quadrant and going
%     counterclockwise.
x = [trv, -trv, -trv, trv];
y = [hrv, hrv, -hrv, -hrv];

%Matrix data for visualisation using 'surf'
Xs = repmat([x, x(1)], [2, 1]);
Ys = repmat([y, y(1)], [2, 1]);
Zs = repmat([0 ; lrv], [1, 5]);
Xp = repmat(x', [1, 2]);
Yp = repmat(y', [1, 2]);
Zp = repmat([0, lrv], [4, 1]);

%Plot
hg(1) = surf(Xs, Ys, Zs);
hg(2) = patch('XData', Xp, 'YData', Yp, 'ZData', Zp);

set(hg, 'FaceColor', clr, 'FaceAlpha', alfa);
end

function hg = drawUserPoint(hAx, hg, x, y, z)
%drawUserPoint Draws a marker at the point described by (x,y,z) in the axes
%'hAx'. If the graphics object 'hg' has already been initialised then
%simply update the 'XData', 'YData', 'ZData' values.

if isempty(hg)
   hg = plot3(hAx, x, y, z, ...
       'Marker'         , 'o', ...
       'MarkerFaceColor', 'b', ...
       'MarkerEdgeColor', 'k', ...
       'Tag', 'BezierControlPoint');
   view(hAx, [0 90]);
else
    set(hg, 'XData', x, 'YData', y, 'ZData', z);
end

end

%Callbacks and UI functionality
function cbMouseClicked(hFig, evt)
%cbMouseClicked Updates the figure app data with the current mouse status
%whenever the user clicks on the figure 'hFig'.

%Covnersion factor from axes units to inches
hg = evt.HitObject;
if ~isa(hg, 'matlab.graphics.axis.Axes')
    hg = ancestor(hg, 'axes');
    if isempty(hg)
        return
    end
end
units = hg.Units;
set(hg, 'Units', 'inches');
w = diff(hg.XLim) / hg.Position(3);
h = diff(hg.YLim) / hg.Position(4);
set(hg, 'Units', units);

%Retrieve the coordinates of the Bezier control points
hCP      = getappdata(hFig, 'BezierControlPointHandle');
cpCoords = [hCP.XData ; hCP.YData ; hCP.ZData]';

%Tolerance for detecting hits on the marker - Based on marker size
sz = hCP.MarkerSize;
tol = mean((sz / 72) .* [w, h]);

%Get coordinates of the mouse click in the axes coordinate system 
xyzIntersection = evt.IntersectionPoint;

%Distance of mouse-click from any of the marker points
r = sqrt(sum((cpCoords - xyzIntersection).^2, 2));

%Only interested in what happends when we click on a control point
idx = (r < tol);
if ~any(idx)
    return
end

%Tell the figure
setappdata(hFig, 'MouseDepressedOnMarker', true);
setappdata(hFig, 'MarkerIndex'           , idx);

%Let the user know we have changed mode
set(hFig, 'pointer', 'hand');

%Next actions depend on what the user does...
%   1. If the user still have the mouse depressed then the
%      'WindowButtonMotionFcn' should fire and the Bezier curves will be
%      updated.
%   2. If the user releases the mouse then nothing will happen.

end

function cbMouseReleased(hFig, evt)
%cbMouseReleased Executes when the user releases the mosue.
%
% Actions:
%   1. Tell the figure appdata that the mouse is no longer on the marker.
%   2. Update the pointer to indicate 'drag-and-drop' mode has ended.

setappdata(hFig, 'MouseDepressedOnMarker', false);
set(hFig, 'pointer', 'arrow');

end

function cbMouseMotion(hFig, evt)
%cbMouseMotion Executes whenever the user moves the mouse on the figure.

%Are we starting from a marker?
bMarker = getappdata(hFig, 'MouseDepressedOnMarker');

if ~bMarker %Escape route
    return
end

%What point are we at?
xyzIntersection = evt.IntersectionPoint;
if any(isnan(xyzIntersection)) %Prevent the user from moving off the axes
    xyzIntersection = getappdata(hFig, 'LastValidIntersection');
end
setappdata(hFig, 'LastValidIntersection', xyzIntersection);

%Get the marker handle and hit index
hCP = getappdata(hFig, 'BezierControlPointHandle');
idx = getappdata(hFig, 'MarkerIndex');

%Update the graphics data
xData = hCP.XData;
yData = hCP.YData;
zData = hCP.ZData;
xData(idx) = xyzIntersection(1);
yData(idx) = xyzIntersection(2);
% zData(idx) = xyzIntersection(3);%Don't update z as we are in the XY plane
set(hCP, 'XData', xData, 'YData', yData); %, 'ZData', zData);

%Calculate the new Bezier curve coordinates
Bi          = getappdata(hFig, 'BernsteinPolynomials');
coords      = [xData ; yData ; zData];
bezierCurve = coords * Bi;

%Apply to the graphics object
hBez = getappdata(hFig, 'BezierCurveHandle');
hCS  = getappdata(hFig, 'CrossSectionPatchHandle');
set(hBez, 'XData', bezierCurve(1, :), 'YData', bezierCurve(2, :), 'ZData', bezierCurve(3, :));
set(hCS , ...
    'XData', [bezierCurve(1, :), bezierCurve(1, 1)], ...
    'YData', [bezierCurve(2, :), bezierCurve(2, 1)], ...
    'ZData', [bezierCurve(3, :), bezierCurve(3, 1)]);

%Update
drawnow

end

%Modified version of built-in 'ginput' function
function [out1,out2,out3] = custom_ginput(arg1)
%GINPUT Graphical input from mouse.
%   [X,Y] = GINPUT(N) gets N points from the current axes and returns
%   the X- and Y-coordinates in length N vectors X and Y.  The cursor
%   can be positioned using a mouse.  Data points are entered by pressing
%   a mouse button or any key on the keyboard except carriage return,
%   which terminates the input before N points are entered.
%
%   [X,Y] = GINPUT gathers an unlimited number of points until the
%   return key is pressed.
%
%   [X,Y,BUTTON] = GINPUT(N) returns a third result, BUTTON, that
%   contains a vector of integers specifying which mouse button was
%   used (1,2,3 from left) or ASCII numbers if a key on the keyboard
%   was used.
%
%   Examples:
%       [x,y] = ginput;
%
%       [x,y] = ginput(5);
%
%       [x, y, button] = ginput(1);
%
%   See also GTEXT, WAITFORBUTTONPRESS.

%   Copyright 1984-2015 The MathWorks, Inc.
%       - Modified for Bezier curve generation by C.Szczyglowski, 16/10/19

out1 = []; out2 = []; out3 = []; y = [];

%EDIT - Initialise point hg-object as empty
hg = [];

if ~matlab.ui.internal.isFigureShowEnabled
    error(message('MATLAB:hg:NoDisplayNoFigureSupport', 'ginput'))
end

% Check Inputs
if nargin == 0
    how_many = -1;
    b = [];
else
    how_many = arg1;
    b = [];
    if  ~isPositiveScalarIntegerNumber(how_many)
        error(message('MATLAB:ginput:NeedPositiveInt'))
    end
    if how_many == 0
        % If input argument is equal to zero points,
        % give a warning and return empty for the outputs.
        warning (message('MATLAB:ginput:InputArgumentZero'));
    end
end

% Get figure
fig = gcf;
drawnow;
figure(gcf);

% Make sure the figure has an axes
hAx = gca(fig);

% Setup the figure to disable interactive modes and activate pointers.
initialState = setupFcn(fig);

% onCleanup object to restore everything to original state in event of
% completion, closing of figure errors or ctrl+c.
c = onCleanup(@() restoreFcn(initialState));

drawnow
char = 0;

while how_many ~= 0
    waserr = 0;
    try
        keydown = wfbp;
    catch %#ok<CTCH>
        waserr = 1;
    end
    if(waserr == 1)
        if(ishghandle(fig))
            cleanup(c);
            error(message('MATLAB:ginput:Interrupted'));
        else
            cleanup(c);
            error(message('MATLAB:ginput:FigureDeletionPause'));
        end
    end
    % g467403 - ginput failed to discern clicks/keypresses on the figure it was
    % registered to operate on and any other open figures whose handle
    % visibility were set to off
    figchildren = allchild(0);
    if ~isempty(figchildren)
        ptr_fig = figchildren(1);
    else
        error(message('MATLAB:ginput:FigureUnavailable'));
    end
    %         old code -> ptr_fig = get(0,'CurrentFigure'); Fails when the
    %         clicked figure has handlevisibility set to callback
    if(ptr_fig == fig)
        if keydown
            char = get(fig, 'CurrentCharacter');
            button = abs(get(fig, 'CurrentCharacter'));
        else
            button = get(fig, 'SelectionType');
            if strcmp(button,'open')
                button = 1;
            elseif strcmp(button,'normal')
                button = 1;
            elseif strcmp(button,'extend')
                button = 2;
            elseif strcmp(button,'alt')
                button = 3;
            else
                error(message('MATLAB:ginput:InvalidSelection'))
            end
        end
        
        if(char == 13) % & how_many ~= 0)
            % if the return key was pressed, char will == 13,
            % and that's our signal to break out of here whether
            % or not we have collected all the requested data
            % points.
            % If this was an early breakout, don't include
            % the <Return> key info in the return arrays.
            % We will no longer count it if it's the last input.
            break;
        end
        
        axes_handle = gca;
        if ~isa(axes_handle,'matlab.graphics.axis.Axes')
            % If gca is not an axes, warn but keep listening for clicks.
            % (There may still be other subplots with valid axes)
            warning(message('MATLAB:Chart:UnsupportedConvenienceFunction', 'ginput', axes_handle.Type));
            continue
        end
        
        drawnow;
        pt = get(axes_handle, 'CurrentPoint');
        how_many = how_many - 1;
                
        out1 = [out1;pt(1,1)]; %#ok<AGROW>
        y = [y;pt(1,2)]; %#ok<AGROW>
        b = [b;button]; %#ok<AGROW>
        
        %EDIT - Draw the point the user has defined
        hg = drawUserPoint(hAx, hg, out1, y, b);
    end
end

% Cleanup and Restore
cleanup(c);

if nargout > 1
    out2 = y;
    if nargout > 2
        out3 = b;
    end
else
    out1 = [out1 y];
end

end

function valid = isPositiveScalarIntegerNumber(how_many)
valid = ~ischar(how_many) && ...            % is numeric
    isscalar(how_many) && ...           % is scalar
    (fix(how_many) == how_many) && ...  % is integer in value
    how_many >= 0;                      % is positive
end

function key = wfbp
%WFBP   Replacement for WAITFORBUTTONPRESS that has no side effects.

fig = gcf;
current_char = []; %#ok<NASGU>

% Now wait for that buttonpress, and check for error conditions
waserr = 0;
try
    h=findall(fig,'Type','uimenu','Accelerator','C');   % Disabling ^C for edit menu so the only ^C is for
    set(h,'Accelerator','');                            % interrupting the function.
    keydown = waitforbuttonpress;
    current_char = double(get(fig,'CurrentCharacter')); % Capturing the character.
    if~isempty(current_char) && (keydown == 1)          % If the character was generated by the
        if(current_char == 3)                           % current keypress AND is ^C, set 'waserr'to 1
            waserr = 1;                                 % so that it errors out.
        end
    end
    
    set(h,'Accelerator','C');                           % Set back the accelerator for edit menu.
catch %#ok<CTCH>
    waserr = 1;
end
drawnow;
if(waserr == 1)
    set(h,'Accelerator','C');                          % Set back the accelerator if it errored out.
    error(message('MATLAB:ginput:Interrupted'));
end

if nargout>0, key = keydown; end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end

function initialState = setupFcn(fig)

% Store Figure Handle.
initialState.figureHandle = fig;

% Suspend figure functions
initialState.uisuspendState = uisuspend(fig);

% Disable Plottools Buttons
initialState.toolbar = findobj(allchild(fig),'flat','Type','uitoolbar');
if ~isempty(initialState.toolbar)
    initialState.ptButtons = [uigettool(initialState.toolbar,'Plottools.PlottoolsOff'), ...
        uigettool(initialState.toolbar,'Plottools.PlottoolsOn')];
    initialState.ptState = get (initialState.ptButtons,'Enable');
    set (initialState.ptButtons,'Enable','off');
end

%Setup empty pointer
cdata = NaN(16,16);
hotspot = [8,8];
set(gcf,'Pointer','custom','PointerShapeCData',cdata,'PointerShapeHotSpot',hotspot)

% Create uicontrols to simulate fullcrosshair pointer.
initialState.CrossHair = createCrossHair(fig);

% Adding this to enable automatic updating of currentpoint on the figure
% This function is also used to update the display of the fullcrosshair
% pointer and make them track the currentpoint.
set(fig,'WindowButtonMotionFcn',@(o,e) dummy()); % Add dummy so that the CurrentPoint is constantly updated
initialState.MouseListener = addlistener(fig,'WindowMouseMotion', @(o,e) updateCrossHair(o,initialState.CrossHair));

% Get the initial Figure Units
initialState.fig_units = get(fig,'Units');
end

function restoreFcn(initialState)
if ishghandle(initialState.figureHandle)
    delete(initialState.CrossHair);
    
    % Figure Units
    set(initialState.figureHandle,'Units',initialState.fig_units);
    
    set(initialState.figureHandle,'WindowButtonMotionFcn','');
    delete(initialState.MouseListener);
    
    % Plottools Icons
    if ~isempty(initialState.toolbar) && ~isempty(initialState.ptButtons)
        set (initialState.ptButtons(1),'Enable',initialState.ptState{1});
        set (initialState.ptButtons(2),'Enable',initialState.ptState{2});
    end
    
    % UISUSPEND
    uirestore(initialState.uisuspendState);
end
end

function updateCrossHair(fig, crossHair)
% update cross hair for figure.
gap = 3; % 3 pixel view port between the crosshairs
cp = hgconvertunits(fig, [fig.CurrentPoint 0 0], fig.Units, 'pixels', fig);
cp = cp(1:2);
figPos = hgconvertunits(fig, fig.Position, fig.Units, 'pixels', fig.Parent);
figWidth = figPos(3);
figHeight = figPos(4);

% Early return if point is outside the figure
if cp(1) < gap || cp(2) < gap || cp(1)>figWidth-gap || cp(2)>figHeight-gap
    return
end

set(crossHair, 'Visible', 'on');
thickness = 1; % 1 Pixel thin lines.
set(crossHair(1), 'Position', [0 cp(2) cp(1)-gap thickness]);
set(crossHair(2), 'Position', [cp(1)+gap cp(2) figWidth-cp(1)-gap thickness]);
set(crossHair(3), 'Position', [cp(1) 0 thickness cp(2)-gap]);
set(crossHair(4), 'Position', [cp(1) cp(2)+gap thickness figHeight-cp(2)-gap]);
end

function crossHair = createCrossHair(fig)
% Create thin uicontrols with black backgrounds to simulate fullcrosshair pointer.
% 1: horizontal left, 2: horizontal right, 3: vertical bottom, 4: vertical top
for k = 1:4
    crossHair(k) = uicontrol(fig, 'Style', 'text', 'Visible', 'off', 'Units', 'pixels', 'BackgroundColor', [0 0 0], 'HandleVisibility', 'off', 'HitTest', 'off'); %#ok<AGROW>
end
end

function cleanup(c)
if isvalid(c)
    delete(c);
end
end

function dummy(~,~)
end
