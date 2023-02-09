classdef (ConstructOnLoad) Component < awi.model.Entity & awi.mixin.FEable & awi.mixin.Mirrorable
    %Component represents part of a model that has position in 3D space and
    %an associated mass. The Component class can be converted to an
    %equivalent Finite Element representation and supports the capability
    %to link objects of the same type and 'mirror' their geometry.
    %
    % The position of a Component can be described in 2 ways:
    %
    %   * Absolute Coordinates - Defining the 'Origin' of the Component
    %   will position the object at the coordinates given by 'Origin' in
    %   the GLOBAL coordinate system. The value of 'Origin' will override
    %   all other position properties.
    %
    %   * Relative Coordinates - This method allows Component classes to be
    %   positioned relative to their parent object.
    %
    %       - The position of the Component with respect to its parent is
    %       controlled by the 'Position' property which, by default, is set
    %       to (0,0,0). That is to say, the component will be coincident
    %       with its parent.
    %
    %       - The properties 'XOffset', 'YOffset' & 'ZOffset' define
    %       offsets in the global coordinate system which allow the
    %       component to be positioned away from the parent.
    %
    %       - If the parent component is also a 'Stick' object then the
    %       property 'SOffset', along with 'SOffsetFlag' can be used to
    %       position the component along the length of the parent 'Stick'.
    %       The axis along which 'SOffset' is defined is controlled by
    %       'SOffsetFlag'. Any offsets in the global coordinate system are
    %       added after interpolation along the parent Stick object.
    %
    % TODO - Rename `SOffset` to `SOffset_eta` to be consistent
    
    properties (AbortSet, SetObservable)
        %Every component has a display name
        DisplayName; %TODO - Get rid
    end
    
    %Describing the position of the Component as a point in 3D space
    properties (AbortSet, SetObservable)
        %Offset in the global X-direction from the parent's 'Position'
        XOffset = 0;
        %Offset in the global Y-direction from the parent's 'Position'
        YOffset = 0;
        %Offset in the global Z-direction from the parent's 'Position'
        ZOffset = 0;
        %Distance along the parent Stick object using the span/axes as
        %defined by 'SOffsetFlag'
        SOffset = 0;
        %Defines which axis or span is to be used when calculating the
        %offset position along the parent Stick object.
        SOffsetFlag  = 'R';
        %Absolute position in 3D space. Overrides all other position
        %properties.
        Origin
    end
            
    %Absolute & Relative Positions
    properties (Dependent)
        %Position of the Component in the parent hgtransform
        Position
        %'AbsPosition' is the absolute position of the origin of a drawable
        %object in the global coordinate system
        AbsPosition
    end
    
    %Appearance of the point - TODO Add set methods 
    properties
        MarkerStyle;
        MarkerSize;
        MarkerFaceColor;
        MarkerEdgeColor
    end
    
    %Mass
    properties (Dependent)
        %This component's mass (aggregated from all PointMass objects associated with this one)
        ThisMass;        
        %Sum of this component's mass, and that of all it's children
        TotalMass;      
        %Centre of Gravity of this object excluding any child objects
        ThisCoG
        %Centre of Gravity of this object including any child objects
        CoG
    end
    
    %Helper properties
    properties (Dependent, Hidden)
        %Helper property for getting around a Matlab bug
        % TODO - Remove this once bug is fixed.
        CoG_
    end
    
    methods % set / get
        function set.XOffset(obj, val)      %set.XOffset
            %set.XOffset Set method for the property 'XOffset'
            %
            %    - 'XOffset' must be a scalar numeric
            
            validateattributes(val, {'numeric'}, {'scalar', 'nonnan', ...
                'finite', 'real'}, 'XOffset', class(obj));
            obj.XOffset = val;
        end
        function set.YOffset(obj, val)      %set.YOffset
            %set.YOffset Set method for the property 'YOffset'
            %
            %    - 'YOffset' must be a scalar numeric
            
            validateattributes(val, {'numeric'}, {'scalar', 'nonnan', ...
                'finite', 'real'}, 'YOffset', class(obj));
            obj.YOffset = val;
        end
        function set.ZOffset(obj, val)      %set.ZOffset
            %set.ZOffset Set method for the property 'ZOffset'
            %
            %    - 'ZOffset' must be a scalar numeric
            
            validateattributes(val, {'numeric'}, {'scalar', 'nonnan', ...
                'finite', 'real'}, 'ZOffset', class(obj));
            obj.ZOffset = val;
        end
        function set.SOffset(obj, val)      %set.SOffset
            %set.SOffset Set method for the property 'SOffset'
            %    - 'obj.SOffset' must be a scalar numeric with a value in the
            %    range [0 : 1]
            validateattributes(val, {'numeric'}, {'scalar', ...
                'nonnegative', '<=', 1}, class(obj), 'Offset');
            obj.SOffset = val;
        end
        function set.SOffsetFlag(obj, val)  %set.SOffsetFlag
            %set.SOffsetFlag Set method for the property 'SOffsetFlag'.
            %
            %   - 'obj.SOffsetFlag' must be one of the following tokens:
            %     {'X', 'Y', 'Z', 'R'}.
            
            validatestring(val, {'X', 'Y', 'Z', 'R'}, class(obj), ...
                'SOffsetFlag');
            obj.SOffsetFlag =  val;
        end
        function set.Origin(obj, val)       %set.Origin
            %set.Origin Set method for the property 'Origin'
            %
            %    - 'Origin' must be a 3 element row vector and all elements
            %      of the vector must be nonnan, finite and real.
            
            validateattributes(val, {'numeric'}, {'row', 'numel', 3, ...
                'nonnan', 'finite', 'real'}, class(obj), 'Origin');
            obj.Origin = val;
        end
        function val = get.Position(obj)    %get.Position
            %get.Position Get method for the dependent property 'Position'
            %
            % 'Position' is the ...
            
            if isempty(obj.Origin) %Origin overrides everything!
                
                %Start from [0, 0, 0]
                val = [0, 0, 0];
                
                if isa(obj.Parent, 'mvc.model.Collector')
                    ph = obj.Parent.Parent;
                else
                    ph = obj.Parent;
                end
                
                if isempty(ph) && isa(obj, 'awi.mixin.BeamPropertyObject')
                    ph = obj.BeamHandle;
                end
                
                %Go along the parent span (if applicable)
                if isa(obj, 'mvc.mixin.Collectable') && isa(ph, 'awi.model.Stick')
                    
                    %Then the offset represents distance along the line of the parent
                    val = val + ph.s2pos(obj.SOffset, obj.SOffsetFlag);
                    
                end
                
                %Add offsets
                val = val + [obj.XOffset, obj.YOffset, obj.ZOffset];
                
            else
                
                %Everything is relative to the parent!
                if obj.HasParent
                    if isa(obj.Parent, 'mvc.model.Collector')
                        phAbsPos = obj.Parent.Parent.AbsPosition;
                    else
                        phAbsPos = obj.Parent.AbsPosition;
                    end
                elseif isa(obj, 'awi.mixin.BeamPropertyObject') && ~isempty(obj.BeamHandle)
                    phAbsPos = obj.BeamHandle.AbsPosition;
                else
                    phAbsPos = [0, 0, 0];
                end
                val = obj.Origin - phAbsPos;
                
            end
        end
        function val = get.AbsPosition(obj) %get.AbsPosition
            %get.AbsPosition Get method for the property 'AbsPosition'.
            %
            % 'AbsPosition' is the absolute position of the origin of a
            % drawable object in the global coordinate system accounting
            % for offsets from the parent object.
            
            %Start from internal value for position
            pos = obj.Position;
            
            %Set parent (ph) & child (ch) objects
            ph = obj.Parent;
            
            %Loop through parent until we reach the top of the tree
            while ~isempty(ph)
                if isa(ph, 'awi.model.Framework') %Top of the hierachy?
                    break
                elseif isa(ph, 'mvc.model.Collector')
                    ph = ph.Parent;
                    continue
                end
                pos(end+1, :) = ph.Position; %#ok<AGROW>
                ph = ph.Parent;
            end
            
            %Sum up all relative positions to find absolute position
            val = sum(pos, 1);
        end
        function val = get.DisplayName(obj) %get.DisplayName
            %get.DisplayName Get method for the dependent property
            %'DisplayName'.
            
            %Get the value
            val = obj.DisplayName;
            
            %If nothing
            if isempty(val)
                
                %Go with name instead
                val = obj.Name;
                
            end
            
        end
        function val = get.ThisMass(obj)    %get.ThisMass
            
            %If we are ourselves a point mass
            if isa(obj, 'awi.model.PointMass')
                
                %Start from here
                val = obj.Mass;
                
            else
                
                %Start from zero
                val = 0;
                
            end
            
            %If we have any children
            if isa(obj, 'mvc.mixin.Collectable') && obj.HasChildren
                
                %Look for any Point Masses
                b = arrayfun(@(x)isa(x, 'awi.model.PointMass'), obj.Children);
                
                %Combine their contributions
                if any(b)
                    val = val + sum(obj.Children(b).Mass);
                end
                
                %It's more complicated than that, because point masses are grouped under a collector,
                % so they are actually grand-children of the object, not children
                if obj.HasGrandChildren
                    
                    %Look for any Point Masses
                    b = arrayfun(@(x)isa(x, 'awi.model.PointMass'), obj.GrandChildren);
                    
                    %Combine their contributions
                    if any(b)
                        val = val + sum([obj.GrandChildren(b).Mass]);
                    end
                    
                end
                
            end
            
        end
        function val = get.TotalMass(obj)   %get.TotalMass
            %get.TotalMass Get method for the dependent property
            %'TotalMass'.
            
            %Start with this mass
            val = obj.ThisMass;
            
            %Add any children
            if isa(obj, 'mvc.mixin.Collectable') && obj.HasChildren
                
                %It would seem reasonable to do this
                %Only applicable to Children that are also Components
                % (because a non-Component doesn't have a Mass property)
                b = arrayfun(@(x)isa(x, 'awi.model.Component'), obj.Children);
                % m = m + sum([obj(b).Children.ThisMass])
                %
                % but MATLAB errors because not happy with recursively accessing a dependent property.
                %
                %So do this instead
                L = flatlist(obj.Children(b));
                if isempty(L)
                    return;
                end
                
                %Find any Point Masses
                b = arrayfun(@(x)isa(x, 'awi.model.PointMass'), L);
                
                %And ask each in flat list directly
                if any(b)
                    val = val + sum([L(b).Mass]);
                end
                
            end
            
        end
        function val = get.ThisCoG(obj)     %get.ThisCoG
            %get.ThisCoG Calculates the Centre-of-Gravity of the component
            %without the influence of any child objects.
            
            %Start with something usable
            val = [0, 0, 0];
            
            if obj.HasChildren %CoG related to point masses
                
                %Get all 'awi.model.PointMass' objects --> Use the
                %collector objects to handle all this!
                b  = arrayfun(@(x) isa(x, 'awi.model.PointMasses'), obj.Children);
                massCollectors = obj.Children(b);
                
                %Escape route
                if isempty(massCollectors)
                    return
                end
                
                %Combine all point masses into a single collection
                pm = vertcat(massCollectors.Children);
                
                %Get the location of the point masses in the global frame
                phAbsPos = obj.AbsPosition; %All point masses have the same parent
                pos = vertcat(pm.Position); %Local
                pos = pos + repmat(phAbsPos, [numel(pm), 1]); %Global
                
                %Vectorise masses
                M       = vertcat(pm.Mass);
                TotMass = sum(M);            %Could use 'ThisMass' but we have already indexed to find the point masses
                M       = repmat(M, [1, 3]); %Avoid explicit array expansion
                
                %Calculate moments about the origin point (0,0,0)
                MomAboutOrigin = sum(pos .* M, 1);
                
                %Calculate the CoG
                val = MomAboutOrigin ./ repmat(TotMass, [1, 3]);
                
            end
            
            %CoG related to non-structural mass (mass x length)
            
            %CoG related to material mass (density x volume)
            
        end
        function val = get.CoG(obj)         %get.CoG
            %get.CoG Calculates the Centre-of-Gravity of this component
            %including the influence of any child objects
            
            %Start with the CoG of just this component
            val = obj.ThisCoG;
            
            if obj.HasChildren
                
                %Grab them and filter out non-compatible objects
                ch = obj.Children;
                b  = arrayfun(@(x) isa(x, 'awi.model.Component'), ch);    %Only 'Components' have mass
                b_ = arrayfun(@(x) ~isa(x, 'awi.model.PointMasses'), ch); %Already accounted for direct child masses!
                
                %Trim the list
                ch = ch(and(b, b_));
                
                if isempty(ch) %Escape route
                    return
                end
                
                %Grab all CoG and mass values for the remaining children
                %and include this object as well
                CG = vertcat(ch.CoG_, obj.ThisCoG);
                M  = vertcat(ch.TotalMass, obj.ThisMass);
                TotMass = sum(M);       %Grab total mass
                M  = repmat(M, [1, 3]); %Avoid implicit expansion for back-compat
                
                %Calculate moments about the origin point (0,0,0)
                MomAboutOrigin = sum((CG .* M), 1);
                
                %Calculate the CoG
                val = MomAboutOrigin ./ repmat(TotMass, [1, 3]);
                
            end
            
        end
        function val = get.CoG_(obj)        %get.CoG_
            %Helper method for circumventing Matlab bug
            
            val = obj.CoG;
            
        end
    end
    
    methods % construction / destruction
        
        function obj = Component(varargin)
            
            %Pass it on
            obj@awi.model.Entity(varargin{:});
            
            %Configure collectables
            obj.addCollectionSpec(@awi.model.PointMass     , 'Point Mass', [], [], [], [], 'Point Masses', false);
            obj.addCollectionSpec(@awi.model.FuelMass      , 'Fuel Mass' , [], [], [], [], 'Fuel Masses' , false);
            obj.addCollectionSpec(@awi.model.CoordSys      , 'Coordinate System');
            obj.addCollectionSpec(@awi.model.BluffBody     , 'Bluff Body');
            obj.addCollectionSpec(@awi.model.LiftingSurface, 'Lifting Surface');
            
            
            %Extend property groups
            obj.addPropertyGroup('Mass and Inertia', ...
                'ThisMass', 'This Mass', 'Mass of this object, not including that of any children', [], ...
                'TotalMass', 'Total Mass', 'Mass of this object, including that of any children', []);
            
            %Extend property groups
            obj.addPropertyGroup('Appearance', ...
                'MarkerStyle'    , 'Style of the marker'      , ...
                'MarkerSize'     , 'Size of the marker'       , ...
                'MarkerFaceColor', 'Colour of the marker face', ...
                'MarkerEdgeColor', 'Colour of the marker edge');
            
            %Add DisplayName to General properties group
            obj.addPropertyGroup('General', ...
                'DisplayName', 'DisplayName');
            
            %Add Position property group
            obj.addPropertyGroup('Position', ...
                'Origin'      , 'Origin'   , 'The origin of the component in the global coordinate system', [], ...
                'XOffset'     , 'X-Offset' , 'Offset from the parent position in the global X-direction'  , [], ...
                'YOffset'     , 'Y-Offset' , 'Offset from the parent position in the global Y-direction'  , [], ...
                'ZOffset'     , 'Z-Offset' , 'Offset from the parent position in the global Z-direction'  , [], ...
                'SOffset'     , 'S-Offset' , 'Normalised distance along the parent ''Stick'' object'      , [], ...
                'SOffsetFlag' , 'S-Offset Flag', 'Defines the axis along which ''S-Offset'' is calculated', []);
            
            %Position only relevant if object is Visible
            obj.setPropertyGroup('Position', 'Visible', @obj.Visible);
            
            %Update the mirrored behaviour
            addMirrorProps(obj, {'XOffset', 'YOffset', 'ZOffset', ...
                'SOffset', 'SOffsetFlag', 'Origin'} , { ...
                @obj.mirrorPosition, @obj.mirrorPosition, ...
                @obj.mirrorPosition, [], [], []});
            obj.setMirrorBehaviour(varargin{:});
            
        end
        
    end
    
    methods % converting to FE model
        function FEModel = convertThisToFE(obj, FEModel, varargin)
            %convertThisToFE Converts the 'awi.fe.Component' to an
            %equivalent 'awi.fe.Node' object and generates the
            %'awi.fe.PointMass' object for any masses in the collection.
            
            %Start with the base class
            FEModel = convertThisToFE@awi.mixin.FEable(obj, FEModel, varargin{:});
            
            %TODO - Work out how to handle spar objects
            if isa(obj, 'awi.model.Spar')
                %No need to do anything further as the superclass methods
                %will handle the conversion for these obejcts
                return
            end
            
            %Check if the object is actually a beam object...
            if isa(obj, 'awi.model.Stick') 
                %No need to generate a node for this component as this will
                %be done at the subclass level by the 'convertThisToFE'
                %method for the 'awi.model.Beam' object.
                n = [];
            else
                %Create a node for this component
                n    = awi.fe.Node;
                n.X  = obj.AbsPosition';
                %Add it to the collection
                addFEData(FEModel, n);
            end
            
            %Make the point mass objects
            if isa(obj, 'mvc.mixin.Collectable') && obj.HasChildren
                
                %It would seem reasonable to do this
                %Only applicable to Children that are also Components
                % (because a non-Component doesn't have a Mass property)
                b = arrayfun(@(x)isa(x, 'awi.model.PointMasses'), obj.Children);                
                if ~any(b)
                    return
                end
                
                %Mass objects are collected in 'awi.model.PointMasses'
                %objects.
                %   - Account for multiple  'awi.model.PointMasses' objects
                %     by using 'vertcat'
                mass = vertcat(obj.Children(b).Children); 

                if isempty(mass)
                    return
                end
                
                %Create the 'awi.fe.PointMass' objects
                m = arrayfun(@(i) awi.fe.PointMass, 1 : numel(mass), ...
                    'Unif', false);
                m = horzcat(m{:});
                
                %Assign a reference to the 'awi.fe.Node' object
                set(m, 'Node', n);
                
                %Set up the mapping betweem the 'awi.fe' and 'awi.model'
                %objects.
                map = {'M', 'I11', 'I21', 'I22', 'I31', 'I32', 'I33' ; ...
                    'Mass', 'Inertia11', 'Inertia12', 'Inertia22', ...
                    'Inertia13', 'Inertia23', 'Inertia33'};
                
                %Assign the data to the 'awi.fe.PointMass' objects
                for i = 1 : size(map, 2)
                    set(m, map(1, i) , num2cell([mass.(map{2, i})], 1)');
                end
                
                %Add it to the collection
                addFEData(FEModel, m);
                
            end
            
        end
    end
    
    methods % visualisation
        
        function hg = drawElement(obj, ht, tag)
            
            ArrowSize = 1;
            
            %Caller supplied tag explicitly ?
            if nargin < 3
                
                %No - apply a tag indicating that we are just drawing some origins
                tag = 'Origins';
                
            end
            
            %Assign owner of this transform (so we can easily find it later)
            set(ht, 'UserData', obj);
            
            %Component has NO particular geometry we can draw,
            % but useful to render as a set of axes,
            % to give visual confirmation of position and orientation
            xd = [0, 1, 0.8, NaN, 1, 0.8] .* ArrowSize;
            yd = [0, 0, 0.1, NaN, 0, -0.1] .* ArrowSize;
            zd = 0 .* xd;
            
            %X-axis in red
            hg{1} = line('Parent', ht, ...
                'XData', xd, ...
                'YData', yd, ...
                'ZData', zd, ...
                'Marker', 'none', ...
                'LineStyle', '-.', ...
                'Color', 'r', ...
                'Tag', tag);
            
            %Y-axis in green
            hg{2} = line('Parent', ht, ...
                'XData', zd, ...
                'YData', xd, ...
                'ZData', yd, ...
                'Marker', 'none', ...
                'LineStyle', '-.', ...
                'Color', 'g', ...
                'Tag', tag);
            
            %Z-axis in blue
            hg{3} = line('Parent', ht, ...
                'XData', yd, ...
                'YData', zd, ...
                'ZData', xd, ...
                'Marker', 'none', ...
                'LineStyle', '-.', ...
                'Color', 'b', ...
                'Tag', tag);
            
        end
        
    end
    
    methods (Access = protected) % mirroring the object
        function mirrorPosition(obj, src, ~)
            %mirrorPosition Listener callback for the 'XOffset', 'YOffset',
            %&'ZOffset' properties of the 'Component'. Will flip the sign
            %of one of the properties based on the current value of
            %'obj.MirrorPlane'.
            
            %What property?
            nam = src.Name;
            
            switch obj.MirrorPlane
                case 'XY'
                    switch nam
                        case {'XOffset', 'YOffset'}
                            obj.(nam) = obj.MirrorParent.(nam);
                        case 'ZOffset'
                            obj.(nam) = -obj.MirrorParent.(nam);
                    end
                case 'XZ'
                    switch nam
                        case {'XOffset', 'ZOffset'}
                            obj.(nam) = obj.MirrorParent.(nam);
                        case 'YOffset'
                            obj.(nam) = -obj.MirrorParent.(nam);
                    end
                case 'YZ'
                    switch nam
                        case {'YOffset', 'ZOffset'}
                            obj.(nam) = obj.MirrorParent.(nam);
                        case 'XOffset'
                            obj.(nam) = -obj.MirrorParent.(nam);
                    end
            end
            
        end
        function mirrorElement(obj)
            %mirrorElement Flips the sign of the (x,y,z) offsets based on
            %the value of the 'MirrorPlane'.
                       
            %Superclass first...
            mirrorElement@awi.mixin.Mirrorable(obj);
            
            %Copy 'SOffset' & the flag
            obj.SOffset     = obj.MirrorParent.SOffset;
            obj.SOffsetFlag = obj.MirrorParent.SOffsetFlag;
            
            %Update the (x,y,z) offsets - Could call the 'mirrorPosition'
            %function several times but that would result in multiple calls
            %to 'switch' statements which would slow down the execution.
            %Better to duplicate a small amount of code... (D.I.E)
            switch obj.MirrorPlane
                case 'XY'
                    obj.XOffset = obj.MirrorParent.XOffset;
                    obj.YOffset = obj.MirrorParent.YOffset;
                    obj.ZOffset = -obj.MirrorParent.ZOffset;
                case 'XZ'
                    obj.XOffset = obj.MirrorParent.XOffset;
                    obj.YOffset = -obj.MirrorParent.YOffset;
                    obj.ZOffset = obj.MirrorParent.ZOffset;
                case 'YZ'
                    obj.XOffset = -obj.MirrorParent.XOffset;
                    obj.YOffset = obj.MirrorParent.YOffset;
                    obj.ZOffset = obj.MirrorParent.ZOffset;
            end
            
        end
    end
    
    methods (Static)
        
        function varargout = interp1(varargin)
            %Custom version of interp1 avoids undesirable errors thrown
            %when interp called during the build process, before all data
            %present and correct

            %First check for empty input arguments 
            if any(cellfun(@isempty, varargin(1:2)))
                
                %Bail out
                [varargout{1:nargout}] = deal([]);
                return;           
                
            end
            
            %Then check for validity of interpolants
            if numel(unique(varargin{1})) == 1
                
                %Bail out
                [varargout{1:nargout}] = deal([]);
                return;               
                
            end
            
            %Consistency of inputs ?
            if numel(varargin{2}) == 1 && numel(varargin{1}) > 1
                varargin{2} = repmat(varargin{2}, size(varargin{1}));
            end

            %Now pass it on
            [varargout{1:nargout}] = interp1(varargin{:});
                      
        end
        
    end
    
    methods % defaults
        function obj = makeDefaultComponent(obj);
        
            
        end
    end
    
end
