classdef (ConstructOnLoad) Drawable < matlab.mixin.SetGet
    %
    % Represents a drawable thing, with a position and orientation in 3-D space.
    %
    % Exactly how it is drawn is not implemented here.
    %
    %
    % TODO - Implement a vectorised drawing method which removes all of the
    % hTransforms. Objective is to minimise number of hg objects.
    
    properties (AbortSet, SetObservable)
                 
        %Useful to be able to specify a thing that is Drawable, but only
        % for the benefit of its children, but does not actually draw itself
        Visible = true;
        
    end
        
    properties (Dependent)
    
        VisibleDescendants;
        
    end
    
    properties (Access = private)
        
        %The transform(s) that holds this drawing
        hTransform;
        
    end
    
    methods % setters / getters
        
        function val = get.VisibleDescendants(obj)
        
            %Start here
            val = obj.Visible;
            
            %Enough ?
            if val
                
                %Yes
                return;
                
            end
            
            %Check children too
            for i = 1:numel(obj.Children)
                
                %Check this child TODO: SHOULD call VisibleDescendants, but mcos treats it as illegal
                val = val || obj.Children(i).Visible;
                
                %Enough ?
                if val
                    
                    %Yes
                    return;
                    
                end

            end
            
        end
        
        function val = get.hTransform(obj)
            
            %Start here
            val = obj.hTransform;
            
            %If nothing
            if isempty(val)
                
                %We're done
                return;
                
            end
            
            %Check validity
            b = arrayfun(@(x)ishghandle(x) && isvalid(handle(x)), val);
            if any(~b)
                
                %Clean it up
                val = val(b);
                
                %And no point persisting invalid content
                obj.hTransform = val;
                
            end
            
        end
        
        function set.Visible(obj, val)

            %Make a note
            obj.Visible = val;

            %Show / hide Geomety property group accordingly
            obj.setPropertyGroup('Geometry', 'Visible', val);

        end

    end
    
    methods % construction / destruction
        
        function obj = Drawable(varargin)
            
            %If properties are managed by Nameable
            if isa(obj, 'mvc.mixin.Nameable')
            
                %Geometry only relevant if object is Visible
                obj.addPropertyGroup('Geometry');
                obj.setPropertyGroup('Geometry', 'Visible', @obj.Visible);
                
                %Add Visible to the General property group
                % (because turning Visible off suppresses display of Geometry group)                
                obj.addPropertyGroup('General', 'Visible', 'Visible');
            
            end
            
        end
        
    end
    
    methods (Sealed) % visualisation
        
        function setProperty(obj, prp, val, par)
            
            %For each element of obj
            for i = 1:numel(obj)
            
                %Get associated transform(s)
                ht = obj(i).hTransform;
                
                %If parent supplied
                if nargin > 2
                
                    %Down-select accordingly
                    b = arrayfun(@(x)ancestor(x, get(par, 'Type')) == par, ht);
                    ht = ht(b);
                    
                end
                
                %For whatever's left
                for j = 1:numel(ht)
                    
                    %Do NOT apply directly to the transform, because doing so includes
                    % any child transforms, which is less than ideal for a hierarchical #
                    % product breakdown when setting 'Selected' or 'Visible'
                    % set(obj.hTransform, prp, val);
                    
                    %Do it this way instead, starting from child objects
                    hc = get(ht(j), 'Children');
                    
                    %Exclude any contained transforms
                    hc(arrayfun(@(x)isa(x, 'matlab.graphics.primitive.Transform'), hc)) = [];
                    
                    %Safety-net for a bug in MATLAB HG
                    b = i_safetyNet(hc, prp, val);
                    
                    %Apply property / value pair
                    set(hc(b), prp, val);
                
                    %For hggroups
                    bg = arrayfun(@(x)strcmp(class(x), 'matlab.graphics.primitive.Group'), hc(b));
                    if any(bg)
                        
                        %Propagate value to children as well
                        set(vertcat(hc(b(bg)).Children), prp, val);
                        
                    end
                    
                end
                
            end
            
            function b = i_safetyNet(hc, prp, val)
            
                %Hope for the best
                b = true(size(hc));
                
                %We're concerned specifically about {'Selected', 'on'}
                if ~strcmpi(prp, 'Selected')
                    return;
                elseif ~strcmpi(val, 'on')
                    return;
                end
                
                % ... applied to Surfaces
                bS = arrayfun(@(x)isa(x, 'matlab.graphics.chart.primitive.Surface'), hc);
                if ~any(bS)
                    return;
                end
                
                % ... whose Data are empty
                bE = cellfun(@isempty, get(hc(bS), {'CData'}));
                if ~any(bE)
                    return;
                end
                
                %Must NOT set Selected on for any surfaces whose Data are empty
                fdx = find(bS);
                b(fdx(bE)) = false;
                
            end
            
        end
        
        function ht = draw(obj, ha, varargin)
        
            %Edit - C.Szczyglowski 28/07/2019 
            %   * Add option to prevent recursing through children
            p = inputParser;
            addParameter(p, 'bRecurse', true, @(x)validateattributes(x, {'logical'}, {'scalar'}));
            parse(p, varargin{:});
            
            %Where ?
            if nargin < 2
                
                %Draw in new axes
                %   - There is a bug in later versions of MATLAB where the
                %   behaviour of uifigures differs from normal figures.
                hf = figure;
                ha = axes('Parent', hf);
                set(hf, 'WindowButtonDownFcn',    @ButtonDownCallback, ...
                  'WindowScrollWheelFcn',   @WindowScrollWheelCallback, ...
                  'KeyPressFcn',            @KeyPressCallback, ...
                  'WindowButtonUpFcn',      @ButtonUpCallback)
      
            validateattributes(ha, {'matlab.graphics.axis.Axes'}, {'scalar'}, class(obj), 'hAx');

                %Make sure have the companion coordinate system!
                hb = obj.makeCoordSystem(ha);             
                
            elseif isa(ha, 'uiextras.Panel')
                
                %Make a note
                hp = ha;
                
                %Check if we already have an axes in this panel
                ha = findobj(hp, 'Type', 'axes');
                
                if numel(ha) > 1
                    idx = ismember({ha.Tag}, '-CoordSys');
                    hb  = ha(idx);
                    ha  = ha(~idx);
                elseif numel(ha) == 1
                    %No extra coordinate system present...make it!
                    hb = obj.makeCoordSystem(ha);  
                end
                    
                %Nothing ?
                if isempty(ha)
                    
                    %Create a new container for all axes - additional 
                    %container required for combination of GLT2 and 3-d 
                    %rotate, legend etc
                    hU = uicontainer('Parent', hp);
                    
                    %Create axes within panel for the main drawing
                    ha = axes('Parent', hU, 'Position', [0 0 1 1]);
                    axis(ha, 'equal');
                    
                    %Make the coordinate system and associated axes
                    hb = obj.makeCoordSystem(ha);  
                    
                    %Add context functionality
                    if isa(obj, 'mvc.mixin.Contextable')
                        context(obj, ha);
                    end
                
                end
                
            elseif isa(ha, 'matlab.graphics.axis.Axes')
                
                %Check if we already have an axes in this panel
                ha = findobj(ha.Parent, 'Type', 'axes');
                
                if numel(ha) > 1
                    idx = ismember({ha.Tag}, '-CoordSys');
                    hb  = ha(idx);
                    ha  = ha(~idx);
                else
                    %No extra coordinate system present...make it!
                    hb = obj.makeCoordSystem(ha);  
                end
                
            end
            
            %Axes passed in ?
            if isa(ha, 'matlab.graphics.axis.Axes')
                
                %Clean sheet (expensive ??)
                cla(ha);
                
                %If we are an Application object at top of tree
                if isa(obj, 'mvc.model.Application') && isa(obj, 'mvc.mixin.Searchable') ...
                        && isempty(obj.findall('isa', 'mvc.mixin.Drawable', '-and', 'Visible', true))
                
                    %Nothing to draw, so create a (hopefully) helpful message
                    str = {obj.Name, ''};
                    if isa(obj, 'mvc.mixin.Serializable')
                        str{end+1} = 'Use File->Open... to load content from a previously saved file';
                        str{end+1} = 'Use File->Recent... to load content from a previous session';
                    end
                    if isa(obj, 'mvc.mixin.Importable')
                        str{end+1} = 'Use File->Import... to import content from file';
                    end
                    if isa(obj, 'mvc.mixin.Collectable')
                        str{end+1} = [newline, 'Select a node, right-click and select "Edit" to interact directly with content'];
                    end
                    if isa(obj, 'mvc.mixin.Documentable')
                        str{end+1} = [newline, 'Use Help->Documents... to access documentation and user guidance'];
                    end
                    
                    %Show us
                    text(0.5, 0.5, str, ...
                        'HorizontalAlignment', 'center', ...
                        'Parent', ha);
                    
                    %Clean it up a bit
                    axis(ha, 'tight');
                    set(ha, 'Visible', 'off');
                    set([hb.Children], 'Visible', 'off');
                    
                else
                    
                    %Pre-label
                    xlabel(ha, 'x');
                    ylabel(ha, 'y');
                    zlabel(ha, 'z');
                    
                    %Show the coordinate system
                    set([hb.Children], 'Visible', 'on');
                    
                    %Link the 'View' of the two axes
                    Link = linkprop([ha, hb], {'View'});
                    setappdata(ancestor(ha, 'matlab.ui.Figure'), 'LinkData', Link);
                   
                    %Add a 'ButtonDownFcn' to allow the user to execute the
                    %default mouse-click behaviour
%                     ha.ButtonDownFcn = @obj.cbOnMouseClick;
                    
                    %Seems sensible to only allow callbacks when the axes is
                    %visible
%                     ha.HitTest       = 'on';
%                     ha.PickableParts = 'visible';  
                    
                end
                
            end

            %Allow for an array of objects to be drawn
            for i = 1:numel(obj)
                
                %Applicable to this object ?
                if ~isa(obj(i), 'mvc.mixin.Drawable')
                    
                    %No
                    continue;
                    
                end
                
                %Draw everything in a transform, assigning owner
                % (so we can easily find it later)
                ht = hgtransform('Parent', ha, 'UserData', obj);
                
                %If this element is flagged as Visible
                if obj(i).Visible
                    
                    %Draw this object
                    drawElement(obj(i), ht);
        
                end
                
                %If we are part of a collection
                if isa(obj(i), 'mvc.mixin.Collectable')
                    
                    %Draw child objects
                    for j = 1:numel(obj(i).Children)
                        
                        %Applicable to this object ?
                        if ~isa(obj(i).Children(j), 'mvc.mixin.Drawable') || ...
                            ~p.Results.bRecurse
                            
                            %No
                            continue;
                            
                        end
                        
                        %Pass it on
                        draw(obj(i).Children(j), ht);
                        
                    end
                    
                end
    
                %Account for objects that don't have a 'Position'
                if isprop(obj(i), 'Position')
                    pos = obj(i).Position;
                else
                    pos = [0, 0, 0];
                end

                %Set position and orientation relative to parent
                M = makehgtform('translate', pos); %, strval{:});
                
                %Apply
                set(ht, 'Matrix', M);
            
                %Make a note of transform
                obj(i).hTransform(end+1) = ht;
                
            end
            
            %Set axis limits to be equal
            if isa(ha, 'matlab.graphics.axis.Axes')
                axis(ha, 'equal');
            end
            
        end
        
    end
    
    methods (Access = private)
        function cbOnMouseClick(obj, src, evt)
            %cbOnMouseClick Callback function linked to the axes
            %'ButtonDownFcn'.
            
            
            
        end
    end
    
    methods % visualisation
        
        function hg = drawElement(~, ~, ~)

            %Return a blank object
            hg{1} = gobjects(1);
            
        end
        
    end
    
    methods (Static)
        
        function hb = makeCoordSystem(ha)
            %makeCoordSystemAxes Creates the axes object and the associated
            %coordinate system in an axes that is a peer of the current
            %axes containing the main drawing.
            
            tag = '-CoordSys';
            
            %Check we are dealing with the correct graphics object
            assert(isa(ha, 'matlab.graphics.axis.Axes'), 'Expected the graphics object to be an axes');
            
            %Grab parent of axes
            hU = ha.Parent;
            
            %Do we already have the coordinate system?
            hb = findobj(hU, 'Tag', tag);
            
            %If not then we need to make it!
            if isempty(hb)
                
                %Create axes within panel for the coordinate system
                hb = axes('Parent', hU, 'Position', [0 0.02 0.2 0.2], ...
                    'Tag', tag);
                
                %Add the coordinate system
                mvc.mixin.Drawable.drawCoordinateSystem(hb);
                
            end
            
        end
        
        function drawCoordinateSystem(ha)
           
            %Parse inputs
            if nargin < 1 || ~isa(ha, 'matlab.graphics.axis.Axes')
                return
            end
            
            %Ensure axes is ready to accept new data
            ha.NextPlot = 'add';
            
            %Define direction vectors
            eX = [0, 0, 0 ; 1, 0, 0];
            eY = [0, 0, 0 ; 0, 1, 0];
            eZ = [0, 0, 0 ; 0, 0, 1];
            
            %Plot the lines and the text
            hL = plot3(ha, ...
                [eX(:, 1), eY(:, 1), eZ(:, 1)], ...
                [eX(:, 2), eY(:, 2), eZ(:, 2)], ...
                [eX(:, 3), eY(:, 3), eZ(:, 3)], ...
                'LineWidth', 2);
            hT = text( ...
                [eX(2, 1), eY(2, 1), eZ(2, 1)], ...
                [eX(2, 2), eY(2, 2), eZ(2, 2)], ...
                [eX(2, 3), eY(2, 3), eZ(2, 3)], ...
                {'X', 'Y', 'Z'}, ...
                'Parent', ha, ...
                'VerticalAlignment', 'bottom');
            
            %Update colours
            set(hL, {'Color'}, {'r' ; 'g' ; 'b'});
            set(hT, {'Color'}, {'r' ; 'g' ; 'b'});
            
            %Update axes appearence
            axis(ha, 'equal');
            set(ha, 'Visible', 'off');
            
        end
        
    end
    
end


function ButtonDownCallback(src, ~)
if strcmp(get(src, 'SelectionType'), 'normal')
% -> the left mouse button is clicked once
% enable the interactive rotation
userData = get(gca, 'UserData');
userData.ppos = get(0, 'PointerLocation');
set(gca, 'UserData', userData)
set(gcf,'WindowButtonMotionFcn',@ButtonMotionCallback)
ButtonMotionCallback(src)   
elseif strcmp(get(src, 'SelectionType'), 'extend')
% -> the left mouse button is clicked once
% enable the interactive rotation
userData = get(gca, 'UserData');
userData.ppos = get(0, 'PointerLocation');
set(gca, 'UserData', userData)
set(gcf,'WindowButtonMotionFcn',@ButtonDragCallback)
ButtonDragCallback(src)
elseif strcmp(get(src, 'SelectionType'), 'open')
% -> the left mouse button is double-clicked
% create a datatip
cursorMode = datacursormode(src);
hDatatip = cursorMode.createDatatip(get(gca, 'Children'));

% move the datatip to the position
ax_ppos = get(gca, 'CurrentPoint');
ax_ppos = ax_ppos([1, 3, 5]);  
% uncomment the next line for Matlab R2014a and earlier
% set(get(hDatatip, 'DataCursor'), 'DataIndex', index, 'TargetPoint', ax_ppos)
set(hDatatip, 'Position', ax_ppos)
cursorMode.updateDataCursors    
end
end
function ButtonMotionCallback(~, ~)
% check if the user data exist
if isempty(get(gca, 'UserData'))
    return
end
% camera rotation
userData = get(gca, 'UserData');
old_ppos = userData.ppos;
new_ppos = get(0, 'PointerLocation');


userData.ppos = new_ppos;
set(gca, 'UserData', userData)

dx = (new_ppos(1) - old_ppos(1))*0.25;
dy = (new_ppos(2) - old_ppos(2))*0.25;
camorbit(gca, -dx, -dy)
end

function ButtonDragCallback(~, ~)
% check if the user data exist
if isempty(get(gca, 'UserData'))
    return
end
% camera rotation
userData = get(gca, 'UserData');
old_ppos = userData.ppos;
new_ppos = get(0, 'PointerLocation');

userData = get(gca, 'UserData');
userData.ppos = new_ppos;
set(gca, 'UserData', userData)


dx = (new_ppos(1) - old_ppos(1))*0.01;
dy = (new_ppos(2) - old_ppos(2))*0.01;
camdolly(gca, -dx, -dy, 0)
end

function WindowScrollWheelCallback(~, eventdata)
% set the zoom facor
if eventdata.VerticalScrollCount < 0
    % increase the magnification
    zoom_factor = 1.05;
else 
    % decrease the magnification
    zoom_factor = 0.95;
end
% camera zoom
camzoom(zoom_factor)
end
function KeyPressCallback(~, eventdata)
% check which key is pressed
if strcmp(eventdata.Key, 'uparrow')
    dx = 0; dy = 0.05;
    camdolly(gca, dx, dy, 0)
elseif strcmp(eventdata.Key, 'downarrow')
    dx = 0; dy = -0.05;
    camdolly(gca, dx, dy, 0)
elseif strcmp(eventdata.Key, 'leftarrow')
    dx = -0.05; dy = 0;
    camdolly(gca, dx, dy, 0)
elseif strcmp(eventdata.Key, 'rightarrow')
    dx = 0.05; dy = 0;
    camdolly(gca, dx, dy, 0)
end

% once again check which key is pressed
if strcmp(eventdata.Key, 'space')
    % restore the original axes and exit the explorer
    userData = get(gcf, 'UserData');
    userData.obj.StopAnimation = true;
end
end
function ButtonUpCallback(~, ~)
% clear the pointer position
    set(gca, 'UserData', [])
end
