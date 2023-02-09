classdef (ConstructOnLoad) LiftingSurface < awi.model.Beam
    %LiftingSurface Defines a generic lifting surface.
    %
    % A lifting surface can have:
    %   - Geometry
    %       + Cross-section (aerofoil coordinates)
    %       +
    %   - Mass & Inertia
    %   - Stiffness
    %   -
    %
    % Parameter Sets
    %   - 'Common' (Common-Set) : Properties common to all other parameter
    %     sets. TODO - Get rid of this!
    %
    %   - 'pSet' (Parameter-Set) : Define the LiftingSurface properties
    %     using the projected angle/span parameterisation method.
    %
    %   - 'cSet' (Coordinate-Set) : Define the LiftingSurface by
    %     explicitly defining the coordinates of the leading and
    %     trailing edge.
    %
    %   - 'sSet' (Sweep-Set) : Special case of the parameter set ('pSet')
    %     where the planform is defined using the LE & TE sweep angles and
    %     the root chord value instead of a set of sweep & chord values.
    %
    %   - 'aSet' (Area-Set) : Define the LiftingSurface using ...
    %
    %   - 'hSet' (Hierachy-Set) : Define the LiftingSurface using the
    %     position of the parent object(s) and offsets.
    %
    % Common properties
    %   - Origin :
    %     Absolute coordinates of the location of the root of the
    %     LiftingSurface stick in 3D space.
    %
    %   - SpanVector :
    %     Axis in the global coordinate system along which the span
    %     direction of the LiftingSurface is defined.
    %
    %   - BeamLoc :
    %     Normalised chordwise location of the beam axis.
    %
    %   - BeamLoc_eta :
    %     Non-dimensional spanwise position of 'BeamLoc'.
    %
    %   - Aerofoil :
    %     Aerodynamic cross-section of the LiftingSurface.
    %
    %   - Aerofoil_eta :
    %     Non-dimensional spanwise position of 'Aerofoil'.
    %
    % Parameter Set
    %   - Sweep :
    %   Angle between the projection of the beam in the global
    %   SpanVector-X plane and the global Y axis
    %
    %   - Dihedral :
    %   Angle between the projection of the beam in the global
    %   SpanVector-Z plane and the global Y axis
    %
    %   - Angle-of-Attack :
    %   Angle between the projection of the beam in the global X-Z
    %   plane and the global X axis
    %
    %   - Sweep Location :
    %   Normalised chordwise position of the sweep angle
    %
    %   - Chord :
    %   Chord length in the global X direction
    %
    % Coordinate Set
    %   - LE :
    %
    %   - TE :
    %
    % Hierachy Set
    %   - TipParent :
    %   Parent object of the tip node.
    %
    %   - TipNormPos :
    %   Normalised distance along the parent span that the tip node is
    %   located at.
    %
    %   - TipFlag :
    %   Defines the span vector which 'TipNormPos' is defined along.
    %
    %   - TipXOffset:
    %   X offset in the global coordinate system of the tip node from the
    %   location defined by 'TipNormPos', 'TipParent' & 'TipFlag'.
    %
    %   - TipYOffset :
    %   Y offset in the global coordinate system of the tip node from the
    %   location defined by 'TipNormPos', 'TipParent' & 'TipFlag'.
    %
    %	- TipZOffset:
    %   Z offset in the global coordinate system of the tip node from the
    %   location defined by 'TipNormPos', 'TipParent' & 'TipFlag'.
    %
    % Sweep Set
    %   - RootChord :
    %   Chord value at the root of the lifting surface.
    %
    %   - LEsweep :
    %   Leading edge sweep angle at the start of each section.
    %
    %   - LEsweep_eta :
    %   Non-dimensionsal spanwise position of 'LEsweep'.
    %
    %   - TEsweep :
    %   Trailing edge sweep angle at the start of each section.
    %
    %   - TEsweep_eta :
    %   Non-dimensionsal spanwise position of 'TEsweep'.
    %
    % TODO - Move the definition of NACA coords into the Cross-Section
    % object.
    % TODO - Add method for generating a wing fold at a specified location
    % along the wing with a specified flare angle
    % TODO - Add method for extruding a wing with a fold (span-extension)
    
    %Appearence
    properties (AbortSet, SetObservable)
        %Edge color of the lines for the 'surf' and 'patch' commands
        EdgeColor = 'none'
        %Face colour for planform
        FaceColor = [0, 1, 1];
        %Transparency for planform
        FaceAlpha = 0.25;
    end
    
    %Common properties
    properties (AbortSet, SetObservable)
        %Axis in the global coordinate system along which the span
        %direction of the LiftingSurface is defined.
        % - Do not assign default value of the 'SpanVector' property as
        % this can cause issues with the Event-Listener behaviour in the
        % 'awi.mixin.Beamable' class. The root cause of this is that the
        % property is 'AbortSet' and this means the 'SpanVector' listener
        % does not fire.
        SpanVector 
        %Aerodynamic cross-section of the LiftingSurface. (Filename or
        %NACA number - e.g. "NACA0012")        
        %Strings are problematic, if backward compat with 2015b to be maintained
        Aerofoil = {'NACA0012', 'NACA0012'};        
        %Non-dimensional spanwise position of 'Aerofoil'
        Aerofoil_eta = [0, 1];
    end
    
    %Coordinate Set
    properties (AbortSet, SetObservable)
        %Location of the leading edge points in the global (x,y,z) system
        LE
        %Location of the trailing edge points in the global (x,y,z) system
        TE
    end
    
    %Parameter Set
    properties (AbortSet, SetObservable)
        %Total length of the beam along the global span vector
        Span = 0;
        %Normalised position of the sweep axis along the chord
        SweepLoc = 0.25;
    end

    %Sweep Set
    properties (AbortSet, SetObservable)
        %Chord value at the root of the lifting surface
        RootChord
    end
    
    %Area Set
    properties (AbortSet, SetObservable)
        %Length of each segment of the lifting surface
        SegLength
        %Area of each segment of the lifting surface
        SegArea
        %Taper of each segment of the lifting surface
        TaperRatio
        %Total projected area of the lifting surface
        SurfaceArea
        %Aspect ratio of the lifting surface
        AspectRatio
    end
    
    %Shadow properties
    properties (Dependent)
        %The finest possible distribution of non-dimensional points along
        %the beam.
        Eta_
        %The spanwise position at all planform property positions
        SpanVec_
        %Normalised chordwise position of the beam line at 'obj.Eta_'.
        BeamLoc_i
        %Coordinates of the edges of the aerodynamic panels.
        PanelCoords = struct();
        %Coordinates of the aerodynamic cross-section at position
        %'AeroEta' ->
        AeroProfile
        %The non-dimensional spanwise position of 'AeroProfile'. It is a
        %combination of 'Eta_' and 'Aerofoil_eta'.
        AeroEta
        %Chord and span-normalised coordinates of the aerofoil profile at
        %position 'Aerofoil_eta'.
        NormalisedProfile
    end
    
    %Mean aerodynamic chord
    properties
        %Value of the Mean aerodynamic chord of the lifting surface 
        MAC
        %(x,y,z) position of the leading edge of the MAC in the local
        %coordinate system
        MACxyzLE
        %(x,y,z) position of the trailing edge of the MAC in the local
        %coordinate system
        MACxyzTE
        %(x,y,z) position of the aerodynamic centre (quarter-chord) of the
        %MAC in the local coordinate system
        MACxyzAC
        %(x,y,z) position of the leading edge of the MAC in the global
        %coordinate system
        MACxyzLE_global
        %(x,y,z) position of the trailing edge of the MAC in the global
        %coordinate system
        MACxyzTE_global
        %(x,y,z) position of the aerodynamic centre (quarter-chord) of the
        %MAC in the global coordinate system
        MACxyzAC_global
    end
    
    %Defining the internal structural layout
    properties
        %Distance between ribs (in the local beam axis system)
        RibPitch      = [];
        %Distance between stringers (in the local beam x-direction)
        StringerPitch = [];
        %Defines the type of mesh for Spar/Rib elements {'uniform',
        %'unstructured'}
        MeshType = 'uniform';
    end

    %wing twist properties
    properties
        Twist = [0,0];
        Twist_eta = [0,1];
    end
 
    %Folding wing tip properties
    properties             
        %Stiffness of the joint in all 3 rotational DOF (Rx, Ry, Rz)
        HingeStiffness   = repmat(1e10, [1, 6]);
        %Define the type of joint {pinned | ball-socket | fully-fixed}
        HingeType = 'pinned';   
        %Logical flag for indicating whether this is a folding wing tip
        IsFoldingWingTip = false;   
        %Flare angle (degrees)
        FlareAngle       = 0;
        %Fold angle (degrees)
        FoldAngle        = 0;
        % fold location
        EtaFold = 0.75;        
    end
    properties (SetAccess = private)
        %Coordinate system of the joint at the fold position
        HingeCoordSytem  = [];
    end
    properties (Dependent)
        %Handle to the child LiftingSurface which is a folding wing tip
        FoldingWingTip
        %Logical flag indicating whether a FWT is a child of this object
        %N.B. Only immediate children are checked, not the whole hierachy
        HasFoldingWingTip
    end
    properties (Constant)
        %Allowable types of hinge for the folding wing tip
        ValidHingeTypes = {'pinned', 'fully-fixed', 'ball-socket'};
    end
    
    %Internal properties
    properties (Hidden = true)
        %X-coordinate of the leading edge at the 'obj.Eta_'
        XLE
        %X-coordinate of the trailing edge at the 'obj.Eta_'
        XTE
        %Number of points along the x-axis for defining the aerofoil
        %coordinates
        NumAerofoilPointsX = 50;
    end
    
    %Shortcuts 
    properties (Dependent)
        %Flag indicating whether the object has any control surfaces
        HasControlSurface
        %Shortcut to all control surfaces, including specific control
        %surface types such as Spoiler, Slat, etc.
        ControlSurfaces
        %Shortcut to all child 'awi.model.Aileron' objects
        Ailerons
        %Shortcut to all child 'awi.model.Flap' objects
        Flaps
        %Shortcut to all child 'awi.model.Slat' objects
        Slats
        %Shortcut to all child 'awi.model.Spoiler' objects
        Spoilers
    end
    
    %Helper properties
    properties (Constant)
        %Key word for all beam properties specific to the lifting surface 
        %that are added to the object during construction.
        BeamPropType = 'Planform Property';
        %Helper property for defining specific control surface types
        %   ( function | name | collection name )
        ControlSurfaceSubTypes = { ...
            @awi.model.Aileron, 'Aileron', 'Ailerons' ; ...
            @awi.model.Flap   , 'Flap'   , 'Flaps'    ; ...
            @awi.model.Slat   , 'Slat'   , 'Slats'    ; ...
            @awi.model.Spoiler, 'Spoiler', 'Spoilers'};
    end
    
    methods % set / get
        function set.SpanVector(obj, val)         %set.SpanVector
            %set.SpanVector Set method for the property 'SpanVector'.
            %
            % 'SpanVector' must be one of the following {'Y', 'Z'}

            val = validatestring(val, {'Y', 'Z'}, class(obj), 'SpanVector');            
            obj.SpanVector = val;
        end
        function set.SweepLoc(obj, val)           %set.SweepLoc
            %set.SweepLoc Set method for the property 'SweepLoc'.
            %
            % Rules :
            %    - 'SweepLoc' must be of type "double" with attributes:
            %      scalar, nonnan, finite, real.
            
            validateattributes(val, {'double'}, {'scalar', 'real', ...
                'finite', 'nonnan', 'nonnegative', '<=', 1}, class(obj), 'SweepLoc');
            obj.SweepLoc = val;
        end
        function set.MeshType(obj, val)           %set.MeshType
            val = validatestring(val, {'uniform', 'unstructured'}, class(obj), 'MeshType');
            obj.MeshType = val;
        end
        function set.RibPitch(obj, val)           %set.RibPitch
           validateattributes(val, {'numeric'}, {'scalar', 'positive'}, ...
               class(obj), 'RibPitch');
           obj.RibPitch = val;
        end
        function set.StringerPitch(obj, val)      %set.StringerPitch
            validateattributes(val, {'numeric'}, {'scalar', 'positive'}, ...
                class(obj), 'StringerPitch');
            obj.StringerPitch = val;
        end
        function val = get.Eta_(obj)              %get.Eta_
            %get.Eta_ Get method for the property 'Eta_'.
            %
            % 'Eta' is the finest distribution required to define the
            % sweep, dihedral and AoA orientation angles as well as the
            % chord.
            
            if isempty(obj.SpanVector)
                flag = 'R';
            else
                flag = obj.SpanVector;
            end
            
            %All 'Planform Properties' should have their '_flag' properties
            %set to be equal to the 'SpanVector' of the 'LiftingSurface'.
            val = getEta(obj, 'Planform Property', 'AxisFlag', flag);
            
        end
        function val = get.SpanVec_(obj)          %get.SpanVec_
            %get.SpanVec_ Get method for the property 'SpanVec_'
            
            if isempty(obj.Span)
                val = [];
                return
            end
            
            val = obj.Eta_ .* obj.Span;
            
        end
        function val = get.BeamLoc_i(obj)         %get.BeamLoc_i
            %get.BeamLoc_i Get method for the property 'BeamLoc_i'
            
            %Objects built using the 'hSet' can only be a straight line so
            %override the calculated value.
            %
            %   TODO - Fix the bug in the method 'calculateBeamLoc' for
            %   objects that are constructed with 'hSet'. Alternatively,
            %   move away from this method of calculating the beam location
            %   and switch to a root location and an angle.
            if strcmpi(obj.ActiveSet, 'hSet')
                val = repmat(obj.BeamLoc(1), size(obj.Eta_));
            else
                val = calculateBeamLoc(obj, obj.Eta_, obj.Chord_, obj.XLE);
            end
        end
        function val = get.PanelCoords(obj)       %get.PanelCoords
            %get.PanelCoords Get method for the property 'PanelCoords'.
            %
            % Each lifting surface can be modelled using a panel method.
            % As these methods rely on linear aerodynamic theory (i.e.
            % small disturbances) the panels must be aligned with the flow
            % direction and cannot be orientated with the mean chord line
            % of the aerofoil. That is to say, they are untwisted.
            %
            % To reduce the burden on the analysis block this object
            % calculates the coordinates of the verticies of the planform
            % which can then be discretized into individual panels.
            %
            % 'PanelCoords' is a structure with three fields ('X', 'Y',
            % 'Z'). A structure is used so that only one 'get' method is
            % required to calculate the coordinates, thereby limiting
            % repetitive calls to the 'LE' and 'TE' properties.
            %
            % The coordinates of the LE and TE are in the global coordinate
            % system. They are NOT relative to the position of the parent
            % and therefore the coordinates  of 'PanelCoords' are also
            % ABSOLUTE coordinates and NOT relative to the parent.
            
            %Preallocate blank structure
            temp = struct('X' , []  , 'Y' , [], 'Z', []);
            val  = struct('LE', temp, 'TE', temp);
            
            %Check if required data is present
            if isempty(obj.LE) || isempty(obj.TE)
                return
            end
            
            %Interpolate LE & TE coordinates to the same spanwise position
            %   - The panels must be streamwise! If we do not interpolate
            %     to the same spanwise position then we will have skewed
            %     panels.
            switch obj.SpanVector
                case 'Y'
                    y = unique([obj.LE(2, :), obj.TE(2, :)], 'stable');
                    xLE = awi.model.Component.interp1(obj.LE(2, :), obj.LE(1, :), y);
                    yLE = y;
                    zLE = awi.model.Component.interp1(obj.LE(2, :), obj.LE(3, :), y);
                    xTE = awi.model.Component.interp1(obj.TE(2, :), obj.TE(1, :), y);
                    yTE = y;
                    zTE = awi.model.Component.interp1(obj.TE(2, :), obj.TE(3, :), y);
                case 'Z'
                    z = unique([obj.LE(3, :), obj.TE(3, :)], 'stable');
                    xLE = awi.model.Component.interp1(obj.LE(3, :), obj.LE(1, :), z);
                    yLE = awi.model.Component.interp1(obj.LE(3, :), obj.LE(2, :), z);
                    zLE = z;
                    xTE = awi.model.Component.interp1(obj.TE(3, :), obj.TE(1, :), z);
                    yTE = awi.model.Component.interp1(obj.TE(3, :), obj.TE(2, :), z);
                    zTE = z;
            end
            
            %Grab planform property values
            aoa  = getBPV(obj, 'AoA');
            bLoc = obj.BeamLoc_i;
            chrd = getBPV(obj, 'Chord');
            
            %Swap NaN for zero 
            aoa(isnan(aoa))   = 0;
            chrd(isnan(chrd)) = 0;
            
            %Calculate offset due to AoA
            zrs = zeros(1, numel(aoa));
            switch obj.SpanVector
                case 'Y'
                    dX_LE = zrs;
                    dX_TE = zrs;
                    dZ_LE = tand(aoa)  .* bLoc       .* chrd;
                    dZ_TE = -tand(aoa) .* (1 - bLoc) .* chrd;
                case 'Z'
                    dX_LE = tand(aoa)  .* bLoc       .* chrd;
                    dX_TE = -tand(aoa) .* (1 - bLoc) .* chrd;
                    dZ_LE = zrs;
                    dZ_TE = zrs;                    
            end

            val.LE.X = xLE - dX_LE;
            val.LE.Y = yLE;
            val.LE.Z = zLE - dZ_LE;            
            val.TE.X = xTE - dX_TE;
            val.TE.Y = yTE;
            val.TE.Z = zTE - dZ_TE;
            
        end
        function val = get.AeroEta(obj)           %get.AeroEta
            %get.AeroEta Get method for the dependent property 'AeroEta'.
            %
            % 'AeroEta' is the non-dimensional spanwise position of the
            % aerofoil cross-sections. It is a vector containing the unique
            % combination of 'Eta_' & 'Aerofoil_eta'.
            
            val = unique([obj.Eta_, obj.Aerofoil_eta]);
        end
        function val = get.NormalisedProfile(obj) %get.NormalisedProfile
            %get.NormalisedProfile Get method for the dependent property
            %'NormalisedProfile'.
            %
            % 'NormalisedProfile' is a structure with three fields
            % ('X', 'Y', 'Z'). A structure is used so that only one 'get'
            % method is required to calculate the coordinates, thereby
            % limiting repetitive calls to any dependent properties.
            %
            % 'X', 'Y' & 'Z' are the normalised coordinates of the aerofoil
            % profile associated with this LiftingSurface object. They are
            % stored in a matrix form in order to enable operations such as
            % interpolation and plotting to be vectorised.
            
            val = [];
            return
            
            %TODO - Fix this
            
            %Pass it on
            val = obj.getNormalisedProfile;
            
            
        end
        function val = get.AeroProfile(obj)       %get.AeroProfile
            %get.AeroProfile Get method for the dependent property
            %'AeroProfile'.
            %
            % 'AeroProfile' is a structure with three fields ('X', 'Y',
            % 'Z'). A structure is used so that only one 'get' method is
            % required to calculate the coordinates, thereby limiting
            % repetitive calls to any dependent properties.
            %
            % 'X', 'Y'& 'Z' are the coordinates of the aerofoil profile
            % associated with this LiftingSurface object. They are stored
            % in a matrix form in order to enable operations such as
            % interpolation and plotting to be vectorised.
            %
            % 'X', 'Y' & 'Z' are relative coordinates. The are positioned
            % in the local coordinate system of the LiftingSurface and are
            % relative to the 'XData', 'YData', 'ZData' properties.
            
            val = [];
            return
            %TODO - Fix this.
            %Pass it on
            val = obj.getAeroProfile;
            
        end
        function val = get.HasControlSurface(obj) %get.HasControlSurface
            %get.HasControlSurface Get method for the property
            %'HasControlSurface';
            
            %Assume yes...
            val = true;
            
            %Simple!
            if isempty(obj.ControlSurfaces)
               val = false; 
            end
            
        end
        function val = get.ControlSurfaces(obj)   %get.ControlSurfaces
            %get.ControlSurfaces Get method for the dependent property
            %'ControlSurfaces'.
            
            %Use the helper function!
            val = grabControlSurface(obj, 'awi.model.ControlSurface');            

        end
        function val = get.Ailerons(obj)          %get.Ailerons
            %get.Ailerons Get method for the dependent property 'Ailerons'.
            
            %Use the helper function
            val = grabControlSurface(obj, 'ailerons');
            
        end
        function val = get.Flaps(obj)             %get.Flaps
            %get.Flaps Get method for the dependent property 'Flaps'.
            
            %Use the helper function
            val = grabControlSurface(obj, 'flaps');
            
        end
        function val = get.Slats(obj)             %get.Slats
            %get.Slats Get method for the dependent property 'Slats'.
            
            %Use the helper function
            val = grabControlSurface(obj, 'slats');
            
        end
        function val = get.Spoilers(obj)          %get.Spoilers
            %get.Spoilers Get method for the dependent property 'Spoilers'.
            
            %Use the helper function
            val = grabControlSurface(obj, 'spoiler');
            
        end
        function val = get.MAC(obj)               %get.MAC
            %get.MAC Get method for the dependent property 'MAC'.
            
            %Pass it on
            [val, ~, ~] = calculateMAC(obj);            
            
        end
        function val = get.MACxyzLE(obj)          %get.MACxyzLE
            %get.MACxyzLE Get method for the dependent property
            %'MACxyzLE'
            
            %Pass it on
            [~, val, ~] = calculateMAC(obj); 
            
        end
        function val = get.MACxyzTE(obj)          %get.MACxyzTE
            %get.MACxyzTE Get method for the dependent property
            %'MACxyzTE'
            
            %Pass it on
            [~, ~, val, ~] = calculateMAC(obj); 
            
        end
        function val = get.MACxyzAC(obj)          %get.MACxyzAC
            %get.MACxyzAC Get method for the dependent property
            %'MACxyzAC'
            
            %Pass it on
            [~, ~, ~, val] = calculateMAC(obj);
            
        end
        function val = get.MACxyzLE_global(obj)   %get.MACxyzLE_global
            %get.MACxyzLE_global Get method for the dependent property
            %'MACxyzLE_global'
            
            %Pass it on
            [~, val, ~] = calculateMAC(obj);
            
            %Account for global coordinate system
            val = val + obj.AbsPosition;
            
        end
        function val = get.MACxyzTE_global(obj)   %get.MACxyzTE_global
            %get.MACxyzTE_global Get method for the dependent property
            %'MACxyzTE_global'
            
            %Pass it on
            [~, ~, val, ~] = calculateMAC(obj);
            
            %Account for global coordinate system
            val = val + obj.AbsPosition;
            
        end
        function val = get.MACxyzAC_global(obj)   %get.MACxyzAC_global
            %get.MACxyzAC_global Get method for the dependent property
            %'MACxyzAC_global'
            
            %Pass it on
            [~, ~, ~, val] = calculateMAC(obj);
            
            %Account for global coordinate system
            val = val + obj.AbsPosition;
            
        end
        function set.HingeStiffness(obj, val)     %set.HingeStiffness
            validateattributes(val, {'numeric'}, {'row', 'numel', 6, ...
                'nonnan', 'finite', 'real'});
            obj.HingeStiffness = val;
        end
        function set.HingeType(obj, val)          %set.HingeType
            val = validatestring(val, obj.ValidHingeTypes, class(obj), 'HingeType');
            obj.HingeType = val;
        end
        function set.IsFoldingWingTip(obj, val)   %set.IsFoldingWingTip
            validateattributes(val, {'logical'}, {'scalar'});
            obj.IsFoldingWingTip = val;
        end
        function set.FlareAngle(obj, val)         %set.FlareAngle
            validateattributes(val, {'numeric'}, {'scalar', 'finite', ...
                'real', 'nonnan', '<', 90}, class(obj), 'FlareAngle');
            obj.FlareAngle = val;
            %Update Child/Parent
            if obj.IsFoldingWingTip && ~isempty(obj.Parent) %#ok<*MCSUP>
                obj.Parent.FlareAngle = val;
            elseif obj.HasFoldingWingTip
                obj.FoldingWingTip.FlareAngle = val;
            end
        end
        function set.FoldAngle(obj, val)          %set.FoldAngle
            validateattributes(val, {'numeric'}, {'scalar', 'finite', ...
                'real', 'nonnan', '<', 90}, class(obj), 'FoldAngle');
            obj.FoldAngle = val;
            %Update Child/Parent
            if obj.IsFoldingWingTip && ~isempty(obj.Parent) %#ok<*MCSUP>
                obj.Parent.FoldAngle = val;
            elseif obj.HasFoldingWingTip
                obj.FoldingWingTip.FoldAngle = val;
            end
        end
        function val = get.FoldingWingTip(obj)    %get.FoldingWingTip
            val = [];
            ch  = obj.Children;
            if isempty(ch)
                return
            end
            LS_ch = find(ch, 'Type', 'LiftingSurface');
            if isempty(LS_ch)
                return
            end
            idx   = [LS_ch.IsFoldingWingTip];
            if any(idx)
                val = LS_ch(idx);
            end
        end
        function val = get.HasFoldingWingTip(obj) %get.HasFoldingWingTip
            val = false;
            if isempty(obj.FoldingWingTip)
                return
            else
                val = true;
            end            
        end
    end
    
    methods % construction / destruction 
        function obj = LiftingSurface(varargin)
            %LiftingSurface Class Constructor performs the following
            %actions:
            %   - Invokes the superclass constructor (awi.model.Component)
            %   - Defines the collection specification
            %   - Updates property groups ('Appearance')
            %   - Defines the parameter sets
            %   - Updates the list of allowable compartments
            
            %Pass it on
            obj@awi.model.Beam(varargin{:});
            
            %Define 'Planform Properties'
            addBeamProperty(obj, 'Chord'   , ...
                'Type', obj.BeamPropType);
            addBeamProperty(obj, 'BeamLoc' , ...
                'Value'    , [0.5, 0.5], ...
                'Variation', 'previous', ...
                'Type'     , obj.BeamPropType);
            addBeamProperty(obj, 'Sweep', ...
                'Variation', 'previous', ...
                'Type'     , obj.BeamPropType);
            addBeamProperty(obj, 'Dihedral', ...
                'Variation', 'previous', ...
                'Type'     , obj.BeamPropType);
            addBeamProperty(obj, 'AoA' , ...
                'Variation', 'linear', ...
                'Type'     , obj.BeamPropType);
            addBeamProperty(obj, 'tc', ...
                'Type', obj.BeamPropType, ...
                'Name', 'Thickness-to-Chord Ratio');  
            addBeamProperty(obj, 'LESweep', ...
                'Type', obj.BeamPropType  , ...
                'Name', 'Leading Edge Sweep');
            addBeamProperty(obj, 'TESweep' , ...
                'Type', obj.BeamPropType   , ...
                'Name', 'Trailing Edge Sweep');
            
            %Configure collectables
            obj.addCollectionSpec(@awi.model.ControlSurface, 'Control Surface', [], [], [], [], 'Generic Control Surfaces', true);
            for iCS = 1 : size(obj.ControlSurfaceSubTypes, 1)
                obj.addCollectionSpec(obj.ControlSurfaceSubTypes{iCS, 1}, ...
                    obj.ControlSurfaceSubTypes{iCS, 2}, [], [], [], [], obj.ControlSurfaceSubTypes{iCS, 3}, true);
            end
            obj.addCollectionSpec(@awi.model.Spar, 'Spar');
            
            %Extend property groups
            obj.addPropertyGroup('Appearance',  ...
                'FaceColor', 'Colour of the patch' , ...
                'FaceAlpha', 'Transparencey of the patch');
            
            %Add Parameter Sets
            obj.addParameterSet('Common', ...
                'DisplayName' , 'Common Set', ...
                'Description' , 'Parameters common to all other methods', ...
                'Precedence'  , inf, ...
                'Parent'      , '' , ... %TODO - Does the parent need to be specified? Surely only for the 'hSet'?
                'SpanVector'  , 'Spanwise axis in the global coordinate system'    , ...
                'Aerofoil'    , 'Aerofoil cross section'                           , ...
                'Aerofoil_eta', 'Non-dimensional spanwise position of ''Aerofoil''', ...
                'BeamLoc'     , 'Normalised chordwise position of the beam'        , ...
                'BeamLoc_eta' , 'Non-dimensional spanwise position of ''BeamLoc''');
            obj.addParameterSet('pSet', ...
                'DisplayName', 'Parametric Set', ...
                'Description', ['Classical wing parameterisation method ', ...
                'based on projected sweep, dihedral and AoA angles']     , ...
                'Precedence'  , 2, ...
                'Span'        , 'Total length of the beam along the global span vector', ... 
                'Sweep'       , 'Sweep angle at the start of each section'         , ...
                'Sweep_eta'   , 'Non-dimensional spanwise position of ''Sweep'''   , ...                
                'SweepLoc'    , 'Normalised chordwise position of the sweep angle' , ...
                'Dihedral'    , 'Dihedral angle at the start of each section'      , ...
                'Dihedral_eta', 'Non-dimensional spanwise position of ''Dihedral''', ...
                'AoA'         , 'Angle-of-Attack at the start of each section'     , ...
                'AoA_eta'     , 'Non-dimensional spanwise position of ''AoA'''     , ...
                'Chord'       , 'Chord length in the global X-direction'           , ...
                'Chord_eta'   , 'Non-dimensional spanwise position of ''Chord''');
            %add chord coordinate system flag -- 'local' vs. 'global'
            obj.addParameterSet('cSet', ...
                'DisplayName' , 'Coordinate Set', ...
                'Precedence'  , 1, ...
                'Description' , ['Provides the absolute coordinates of '    , ...
                'the leading edge and trailing edge of the lifting surface'], ...
                'LE', 'Leading edge coordinates (Absolute)', ...
                'TE', 'Trailing edge coordinates (Absolute)');
            obj.addParameterSet('hSet', ....
                'DisplayName', 'Hierachy Set', ...
                'Precedence' , 3, ...
                'Description', ['Defines a straight LiftingSurface ', ...
                'object using an offset along a given parent object.'], ...
                'Parent'      , 'Parent object of the root node', ...
                'TipParent'   , 'Parent object of the tip node', ...
                'TipNormPos'  , 'Normalised location of the tip node along the parent object', ...
                'TipFlag'     , 'The span vector which ''TipNormPos'' is defined along', ...
                'TipXOffset'  , ['X offset in the global coordinate system '   , ...
                'of the root node from the location defined by ''TipNormPos''' , ...
                ', ''TipParent'' & ''TipFlag'''], ...
                'TipYOffset'  , ['Y offset in the global coordinate system '   , ...
                'of the root node from the location defined by ''TipNormPos''' , ...
                ', ''TipParent'' & ''TipFlag'''], ...
                'TipZOffset'  , ['Z offset in the global coordinate system '   , ...
                'of the root node from the location defined by ''TipNormPos''' , ...
                ', ''TipParent'' & ''TipFlag'''], ...
                'Chord'       , 'Chord length in the global X-direction'           , ...
                'Chord_eta'   , 'Non-dimensional spanwise position of ''Chord''');
            obj.addParameterSet('sSet', ...
                'DisplayName', 'Sweep Set', ...
                'Precedence' , 5, ...
                'Description', ['Defines the planform using the leading ', ...
                'and trailing edge sweep angles and the root chord.'], ...
                'Span'        , 'Total length of the beam along the global span vector', ...
                'RootChord'   , 'Chord value at the root of the wing', ...
                'LESweep'     , ['Leading edge sweep angle at the '  , ...
                'start of each section'], ...
                'LESweep_eta' , 'Non-dimensional spanwise position of ''LEsweep''', ...
                'TESweep'     , ['Trailing edge sweep angle at the ', ...
                'start of each section'] , ...
                'TESweep_eta' , 'Non-dimensional spanwise position of ''TEsweep''');
            obj.addParameterSet('aSet', ...
                'DisplayName', 'Area Set', ...
                'Precedence' , 4, ...
                'Description', ['Defines the planform using the aspect ', ...
                'ratio , span and taper ratio.'], ...
                'Span'       , 'Total length of the beam along the global span vector', ... 
                'SegLength'  , 'Length of each segment of the lifting surface', ...
                'AspectRatio', 'Aspect Ratio (AR) of the lifting surface'     , ...
                'TaperRatio' , 'Taper of each segment of the lifting surface' , ...
                'SurfaceArea', 'Total projected area of the lifting surface'  , ...
                'Chord'      , 'Chord length in the global X-direction'       , ...
                'Chord_eta'  , 'Non-dimensional spanwise position of ''Chord''', ...
                'RootChord'  , 'Chord value at the root of the wing');
            
%             %Update the mirrored behaviour
%             addMirrorProps(obj, {'Chord', 'Chord_eta', 'Sweep', ...
%                 'Sweep_eta', 'AoA', 'AoA_eta', 'Dihedral'     , ...
%                 'Dihedral_eta', 'SweepLoc', 'Span'} , {[], [], [], [], ... 
%                 [], [], [], [], [], @obj.assignNegativeMirrorValue});
%             obj.setMirrorBehaviour(varargin{:});
            
            %Update compartments
            updateCompartmentType(obj, 'awi.model.FuelTank');
            
            %Add default builds
            addDefaultBuild(obj, 'HodgesWing', @awi.model.LiftingSurface.makeHodgesWing, '16m semi-span, straight wing');
            addDefaultBuild(obj, 'A320Wing'  , @awi.model.LiftingSurface.makeA320Wing  , 'Simple wing based on Airbus A320');
        end        
    end
    
    methods % visualisation 
        function hg = drawElement(obj, ht)
            
            %Start with base-class
            hg = drawElement@awi.model.Beam(obj, ht);
            
            %Position of origin relative to parent
            O = obj.AbsPosition;
            
            %Plot aerodynamic panel outline
            panelCoords = obj.PanelCoords;  % use get method once!
            xAP = [panelCoords.LE.X ; panelCoords.TE.X] - O(1);
            yAP = [panelCoords.LE.Y ; panelCoords.TE.Y] - O(2);
            zAP = [panelCoords.LE.Z ; panelCoords.TE.Z] - O(3);
            hg{end+1} = surf(xAP, yAP, zAP, 'Parent', ht, ...
                'FaceColor', obj.FaceColor, ...
                'FaceAlpha', obj.FaceAlpha, ...
                'Tag'      , 'Aero Panel');
            
            %Plot the MAC
            %   - Annotate the MAC chord on the lifting surface
            %   - Discrete marker at the MAC quarter-chord
            [~, MACxyzLE_, MACxyzTE_, MACxyzAC_] = calculateMAC(obj);
            clr = [1 0 0];
            hg{end + 1} = line( ...
                'XData' , [MACxyzLE_(1) ; MACxyzTE_(1)], ...
                'YData' , [MACxyzLE_(2) ; MACxyzTE_(2)], ...
                'ZData' , [MACxyzLE_(3) ; MACxyzTE_(3)], ...
                'Parent', ht    , ...
                'Tag'   , 'MAC' , ...
                'LineStyle', '-', ...
                'LineWidth', 2  , ...
                'Color'    , clr);
            hg{end + 1} = line( ...
                'XData', MACxyzAC_(1)    , ...
                'YData', MACxyzAC_(2)    , ...
                'ZData', MACxyzAC_(3)    , ...
                'LineStyle'      , 'none', ...
                'Marker'         , 'o'   , ...
                'MarkerFaceColor', clr   , ...
                'MarkerEdgeColor', 'k'   , ...
                'Parent'         , ht    , ...
                'Tag'            , 'MAC');
            
        end        
    end
    
    methods % class building
        %Note on class building for the awi.model.LiftingSurface class.
        %
        %   * Each method must ensure that the following properties have
        %   been defined:
        %
                
        function build_pSet(obj)
            %build_pSet Builds the LiftingSurface object using the
            %variables belonging to the 'pSet'.
                        
            %Invoke get methods once
            chrd = obj.Chord_;
            sVec = obj.SpanVec_;
            
            %Calculate coordinates of ref. line based on each angle and eta
            sCoords = obj.angle2Beam(abs(sVec), obj.Sweep_);
            dCoords = obj.angle2Beam(abs(sVec), obj.Dihedral_);
            
            %Beam-axis/Sweep-axis offset distance @ root
            LERootOffset = abs(obj.SweepLoc - obj.BeamLoc(1)) .* chrd(1);
            
            %Currently the reference line is located at 'SweepLoc'. Need to
            %find the coordinates of the leading edge in the zero-AoA plane
            %and then work backwards to calculate the coordinates of the
            %beam line.
            obj.XLE = sCoords - (obj.SweepLoc .* chrd) - LERootOffset;
            obj.XTE = sCoords + ((1-obj.SweepLoc) .* chrd) - LERootOffset;
            
            %Adjust the position of the beamline based on the normalised 
            %sweep location
            sCoords = sCoords + (obj.BeamLoc_i - obj.SweepLoc) .* chrd - LERootOffset;
            
            %Assign 'Lineable' properties (XData, YData, ZData)
            switch obj.SpanVector
                case 'X'
                    obj.XData = sVec;
                    obj.YData = sCoords;
                    obj.ZData = dCoords;
                case 'Y'
                    obj.XData = sCoords;
                    obj.YData = sVec;
                    obj.ZData = dCoords;
                case 'Z'
                    obj.XData = sCoords;
                    obj.YData = dCoords;
                    obj.ZData = sVec;
            end
            
        end
        
        function build_cSet(~)
            %build_cSet Build method for the parameter set 'cSet'.
            %
            % This method does not actually do anything as the cSet
            % parameters ('LE' & 'TE') are already calculated by default in
            % the method 'build_pSet'.
            
        end
        
        function build_aSet(~)
            %build_cSet Build method for the parameter set 'aSet'.
            
        end
        
        function build_hSet(obj)
            %build_hSet Builds the LiftingSurface object using the
            %variables belonging to the 'hSet'.
            
            %Invoke superclass method from awi.model.Connector to calculate
            %root and tip coordinates.
            getConnectorXYZ(obj);
            
            %Interpolate the data so that the no. of points match 'Eta_'
            %   - Assume linear distribution
            eta = obj.Eta_; % invoke get method once!
            obj.XData = awi.model.Component.interp1([0, 1], obj.XData, eta);
            obj.YData = awi.model.Component.interp1([0, 1], obj.YData, eta);
            obj.ZData = awi.model.Component.interp1([0, 1], obj.ZData, eta);
            
            %Translate coordinates so that XData, YData & ZData always start
            %at (0,0,0);
            interpPos = s2pos(obj.Parent, obj.SOffset, obj.SOffsetFlag);
            obj.XData = obj.XData - obj.XOffset - interpPos(1);
            obj.YData = obj.YData - obj.YOffset - interpPos(2);
            obj.ZData = obj.ZData - obj.ZOffset - interpPos(3);
            
            %Assign the property 'obj.XLE'.
            obj.XLE = obj.XData - (obj.BeamLoc_i .* obj.Chord_);
        end
        
        function build_sSet(obj)
            %build_sSet Build method for the parameter set 'sSet'.
                        
            %Invoke get method once
            sVec = obj.SpanVec_;
            
            %Calculate the variation in the global X-direction of the
            %leading and trailing edge due to sweep.
            leCoords = obj.angle2Beam(sVec, obj.LESweep_);
            teCoords = obj.angle2Beam(sVec, obj.TESweep_);
            
            %Calculate the variation in the global Z-direction due to
            %dihedral
            dCoords = obj.angle2Beam(abs(sVec), obj.Dihedral_);
            
            %Position LE & TE in global X-direction
            obj.XLE = leCoords - obj.BeamLoc(1) .* obj.RootChord;
            obj.XTE = teCoords + (1 - obj.BeamLoc(1)) .* obj.RootChord;
            
            %Calculate the chord distribution
            obj.Chord = abs(obj.XTE - obj.XLE);
            obj.Chord_eta = obj.Eta_;
            
            %Calculate the position of the beam reference line in global X
            sCoords = leCoords + (obj.BeamLoc_i .* obj.Chord) - (obj.BeamLoc_i(1) * obj.RootChord);
            
            % assign 'Lineable' properties (XData, YData, ZData)
            switch obj.SpanVector
                case 'X'
                    obj.XData = obj.SpanVec_;
                    obj.YData = sCoords;
                    obj.ZData = dCoords;
                case 'Y'
                    obj.XData = sCoords;
                    obj.YData = obj.SpanVec_;
                    obj.ZData = dCoords;
                case 'Z'
                    obj.XData = sCoords;
                    obj.YData = dCoords;
                    obj.ZData = obj.SpanVec_;
            end
                        
        end
        
        function update_pSet(obj)
            %update_pSet Derives the properties belonging to the 'pSet'
            %using the available data in 'obj'.
            %
            % Properties updated:
            %   - Span
            %   - Sweep
            %   - Sweep_eta
            %   - Dihedral
            %   - Dihedral_eta
            %   - AoA
            %   - AoA_eta
            
            %Grab coordinates of stick
            coords = [obj.XData ; obj.YData ; obj.ZData];
            
            if isempty(coords) %Escape route
                return
            end
            
            %How many spanwise positions are we considering?
            nPos = numel(obj.XData);
            
            %Change in coordinates
            dx  = abs(diff(obj.XData));
            dy  = abs(diff(obj.YData));
            dz  = abs(diff(obj.ZData));
                    
            %Calculate the 'Span', 'Sweep', 'Dihedral' & 'AoA' properties
            %   - 'Span' is (tip - root) in span axes
            %   - 'Sweep' is the angle between the projection of the beam
            %   in the global SpanVector-X plane and the global Y axis
            %   - 'Dihedral' is the angle between the projection of the
            %   beam in the global SpanVector-Z plane and the global Y axis
            %   - 'AoA' is the angle between the projection of the beam in
            %   the global X-Z plane and the global X axis
            switch obj.SpanVector
                case 'X'
                    %The 'pSet' angles don't really make sense when we
                    %consider the 'X' spanvector - set to zero for now.
                    obj.Span = coords(1, end) - coords(1, 1);
                    eta      = obj.XData ./ obj.Span;
                    zrs      = zeros(1, nPos); 
                    swp = zrs;
                    dih = zrs;
                    aoa = zrs;
                case 'Y'
                    obj.Span = coords(2, end) - coords(2, 1);
                    eta      = obj.YData ./ obj.Span;                    
                    swp = atand(dx ./ dy);
                    dih = atand(dz ./ dy);
                    aoa = zeros(size(dih)); 
                    %Calculate chord (Get X-position of LE and TE)
                    y = unique([abs(obj.LE(2, :)), abs(obj.TE(2, :))]);
                    xLE = interp1(abs(obj.LE(2, :)), obj.LE(1, :), y);
                    xTE = interp1(abs(obj.TE(2, :)), obj.TE(1, :), y);
                    chord = xTE - xLE;
                    chord_eta = (coords(2, :) - coords(2, 1)) ./ obj.Span;
                case 'Z'
                    obj.Span = coords(3, end) - coords(3, 1);
                    eta      = obj.ZData ./ obj.Span;
                    zrs      = zeros(1, nPos); 
                    swp = atand(dx ./ dz);
                    dih = zrs;
                    aoa = zrs;
                    %Calculate chord (Get X-position of LE and TE)
                    z = unique([obj.LE(3, :), obj.TE(3, :)]);
                    xLE = interp1(obj.LE(3, :), obj.LE(1, :), z);
                    xTE = interp1(obj.TE(3, :), obj.TE(1, :), z);
                    chord = xTE - xLE;
                    chord_eta = (coords(3, :) - coords(3, 1)) ./ obj.Span;
            end
            
            %Assign eta distribution
            obj.Sweep_eta    = eta;
            obj.Dihedral_eta = eta;
            obj.AoA_eta      = eta;
            obj.Chord_eta    = chord_eta;
            
            %Assign dat - Make sure the number of angles matches the number
            %of eta positions
            obj.Chord    = chord;
            obj.Sweep    = [swp, swp(end)];
            obj.Dihedral = [dih, dih(end)];
            obj.AoA      = [aoa, aoa(end)];
                    
        end
        
        function update_cSet(obj)
            %update_cSet Derives the properties belonging to the 'cSet'
            %using the available data in 'obj'.
            %
            % Properties updated:
            %   - LE : Absolute coordinates of the LE in the global
            %          coordinate system.
            %   - TE : Absolute coordinates of the TE in the global
            %          coordinate system.
            
            %Grab origin of the the compnent
            Origin = obj.AbsPosition;
            
            %Grab planform parameters
            aoa  = obj.AoA_;
            bLoc = obj.BeamLoc_i;
            chrd = obj.Chord_;
                                    
            %Calculate z-offset due to twist for LE & TE
            zrs = zeros(1, numel(aoa));
            switch obj.SpanVector
                case 'Y'
                    dX_LE = zrs;
                    dX_TE = zrs;
                    dZ_LE = tand(aoa)  .* bLoc       .* chrd;
                    dZ_TE = -tand(aoa) .* (1 - bLoc) .* chrd;
                case 'Z'
                    dX_LE = tand(aoa)  .* bLoc       .* chrd;
                    dX_TE = -tand(aoa) .* (1 - bLoc) .* chrd;
                    dZ_LE = zrs;
                    dZ_TE = zrs;
            end
            
            %Define leading & trailing edge coordinates (absolute)
            obj.LE = [ ...
                (obj.XData - (bLoc .*chrd)) + Origin(1) + dX_LE; ...
                obj.YData + Origin(2) ; ...
                obj.ZData + Origin(3) + dZ_LE];
            obj.TE = [ ...
                (obj.XData + ((1-bLoc) .* chrd)) + Origin(1) + dX_TE; ...
                obj.YData + Origin(2) ; ...
                obj.ZData + Origin(3) + dZ_TE];
        end
        
        function update_aSet(obj)
            %update_aSet Derives the properties belonging to the 'aSet'
            %using the available data in 'obj'.
            %
            % Properties updated:
            %   - 'SegLength' 
            %   - 'TaperRatio'
            %   - 'AspectRatio'
            %   - 'SurfaceArea'
                        
            %Invoke get methods once!
            sVec = obj.SpanVec_;
            chrd = obj.Chord_;
                        
            %Need chord at start and end of each section
            c0 = chrd(1 : end - 1);
            cN = chrd(2 : end);
            
            %Taper ratio is the ratio of chord values at the start and end
            %of each segment
            obj.TaperRatio = cN ./ c0;
            
            %Segment length is just the length between kinks
            obj.SegLength = diff(abs(sVec));
            
            %Calculate area of each segment
            %   - Use area of a trapezoid
            %   - TODO : Update this so that it uses the shoelace formula            
            obj.SegArea = 0.5 .* (c0 + cN) .* obj.SegLength;
                                    
            %Total surface area is just the sum of the individual segments
            obj.SurfaceArea = sum(obj.SegArea); 
            
            %Aspect ratio is trivial
            %   - Multiply by 2 because we are assuming the other lifting
            %   surface is a mirror of this one
            obj.AspectRatio = 2.* (obj.Span .^2 ./ obj.SurfaceArea);
            
        end
        
        function update_hSet(~)
            %update_hSet Derives the properties belonging to the 'hSet'
            %using the available data in 'obj'.
            %
            % Properties updated:
            %   - NA
            %
            % * *    No properties are updated in this method as the   * *
            % * *    'hSet' is extremely specific in its application   * *
            
        end
        
        function update_sSet(obj)
            %update_sSet Derives the properties belonging to the 'sSet'
            %using the available data in 'obj'.
            %
            % Updates the following parameters:
            %   - RootChord
            %   - LESweep
            %   - LESweep_eta
            %   - TESweep
            %   - TESweep_eta
            
            %Root chord is simple
            obj.RootChord = obj.Chord(1);
            
            %Invoke get method once!
            xLE = obj.XLE;
            eta = obj.Eta_;
            
            %Calculate the leading edge and trailing edge sweep            
            dxLE = diff(xLE); 
            xTE  = xLE + obj.Chord_;
            dxTE = diff(xTE); 
            ds   = obj.SegLength; %Change in span (positive)
            leSweep = atand(dxLE ./ ds);
            teSweep = atand(dxTE ./ ds);                        
            
            %LE/TE sweep is defined at the same eta positions as the beam
            %sweep
            obj.LESweep_eta = obj.Sweep_eta;
            obj.TESweep_eta = obj.Sweep_eta;                        
            
            %Check if we have a single value of sweep
            if numel(leSweep) < numel(eta)
                %Append the sweep vectors with the end value to signify the
                %sweep value at eta = 1;
                leSweep = [leSweep, leSweep(end)];
                teSweep = [teSweep, teSweep(end)];
            end
            
            %'leSweep' and 'teSweep' are currently defined at 'obj.Eta_'
            %because it is based on 'obj.XLE' - Need to interpolate to the
            %sweep eta positions
            obj.LESweep = interp1(eta, leSweep, obj.Sweep_eta);
            obj.TESweep = interp1(eta, teSweep, obj.Sweep_eta);
            
            %TODO - The dynamic property listeners do not seem to be firing
            
        end
                
    end
        
    methods (Access = protected) % converting to FE
        %1D beam model
        function FEModel = generate1DBeamModel(obj, FEModel, varargin)
            %generate1DBeamModel Generates the 'awi.fe.' data necessary to
            %describe a 1D beam-stick representation of a LiftingSurface
            %object.
            
            %Points from geometry
            xG   = obj.XData;
            yG   = obj.YData;
            zG   = obj.ZData;
            %rG   = obj.RData;
            etaR = getEta(obj, 'Planform Property', 'AxisFlag', obj.SpanVector);

            %Pass to 'awi.model.Beam' method...
            FEModel = generate1DBeamModel@awi.model.Beam(obj, FEModel, ...
                'etaR', etaR, 'xR', xG, 'yR', yG, 'zR', zG);            
            
            %Don't bother proceeding if we don't need to model the
            %aerodynamic surface but proceed if the user is forcing the
            %planform to be generated.
            if isempty(FEModel.Nodes)
                return
            end
            if ~obj.ModelAero && ~obj.GenerateLiftSurfPlanform
                return
            end
            
            %Add FE objects to describe the planform
            %   - At this point, the FE model contains data belonging to
            %   the beam.
            %   - Assume that all nodes lie along the beam and we want to
            %   create additional nodes at the LE & TE for the spline.
            %   - If there are additional nodes that do not lie along the
            %   beam then this section will fail as the indexing will not
            %   be correct and additional RBE nodes will be created!
                         
            %Check for control surfaces
            ControlSurf = grabControlSurface(obj, 'awi.model.ControlSurface');
            ControlSurf = checkControlSurfaceOverlap(obj, ControlSurf);
            
            %Calculate the panel coordinates once! (Remember this is in the
            %global frame)
            %   - TODO : What if this is empty?
            PanelCoord = obj.PanelCoords;
            %% Make the LE & TE nodes & Rigid Connections
            
            %How many grid points have already been defined for this FE
            %component? 
            BeamNodes = FEModel.BeamNodes;
            BeamNodes = [BeamNodes(1, :), BeamNodes(2, end)];
            beamX     = [BeamNodes.X];
            nGrid     = size(beamX, 2);
            
            %Preallocate
            NodesLE = arrayfun(@(~) awi.fe.Node    , 1 : nGrid   , 'Unif', false);
            NodesTE = arrayfun(@(~) awi.fe.Node    , 1 : nGrid   , 'Unif', false);
            RB      = arrayfun(@(~) awi.fe.RigidBar, 1 : nGrid   , 'Unif', false);
            %Set     = awi.fe.StructuralSet;                       
            NodesLE = horzcat(NodesLE{:});
            NodesTE = horzcat(NodesTE{:});
            RB      = horzcat(RB{:});
            %
            CoordLE = zeros(3, nGrid);
            CoordTE = zeros(3, nGrid);
            
            %LE & TE coordinates of the lifting surface
            LE_ = [PanelCoord.LE.X ; PanelCoord.LE.Y ; PanelCoord.LE.Z];
            TE_ = [PanelCoord.TE.X ; PanelCoord.TE.Y ; PanelCoord.TE.Z];
            
            %Define coordinates of LE & TE at the beam coordinate positions
            switch obj.SpanVector
                case 'X'
                    error('Update code for ''SpanVector'' = X');
                case 'Y'
                    ind = 2;
                case 'Z'
                    ind = 3;
            end
            CoordLE(1, :)      = interp1(LE_(ind, :), LE_(1, :), beamX(ind, :));
            CoordTE(1, :)      = interp1(TE_(ind, :), TE_(1, :), beamX(ind, :));  
            CoordLE([2, 3], :) = [beamX(2, :) ; beamX(3, :)];
            CoordTE([2, 3], :) = [beamX(2, :) ; beamX(3, :)];
            
            %Populate 'awi.model.fe.Nodes' object(s)
            set(NodesLE, {'X'}, num2cell(CoordLE, 1)');
            set(NodesTE, {'X'}, num2cell(CoordTE, 1)');
            
            %Define the rigid bar independant/dependant nodes
            set(RB, {'NodesI'}, num2cell(BeamNodes'));
            set(RB, {'NodesD'}, num2cell([NodesLE ; NodesTE], 1)');
            
            %Define dependent DOFs for the rigid bar(s)
            set(RB, 'CN', 123456);
            
            %Put all the nodes into a single collection
            %Set.Nodes = horzcat(FEModel.Nodes, NodesLE, NodesTE)';    
            
            %Add objects to the FE model
            % - Add the aero control surfaces seperately to other panels
            addFEData(FEModel, NodesLE, NodesTE, RB); %, Set);
            
            %% Define aerodynamic panels
            
            %No need to proceed if we aren't modelling the aero
            if ~obj.ModelAero
                return
            end
            
            %Create a single 'awi.fe.AeroProp' object
            AeroProp = awi.fe.AeroProp;  
            
            %Define the aerodynamic panels
            if isempty(ControlSurf) || ~obj.ModelControlSurf
                AeroPanel     = obj.defineAeroPanelsNoControlSurf(PanelCoord);  
                ControlPanels = [];
                ControlSurfs  = [];
            else
                %Can only mesh control surfaces that have a 'Label'
                assert(~any(cellfun(@isempty, {ControlSurf.Label})) , ...
                    ['Each control surface must have its ''Label'' ', ...
                    'property defined in order to convert it to a ' , ...
                    'valid MSC.Nastran model.']);
                %Define the panels belonging to just the lifting surface
                AeroPanel     = obj.defineAeroPanelsWithControlSurf(ControlSurf, PanelCoord);
                %Define the panels & other FE objects belonging to the
                %control surfaces
                if obj.ModelControlSurfStructure
                    %Define the full-suite of FE objects including
                    %structural and aerodynamic objects
                    arrayfun(@(cs) convertThisToFE(cs, FEModel, BeamNodes), ControlSurf, 'Unif', false);
                    ControlPanels = [];
                else
                    %Just define the aerodynamic objects and a coordinate
                    %system
                    ControlPanels = obj.defineControlSurfacePanels(ControlSurf);
                    hingeCoordSys = arrayfun(@(cs) defineHingeCoordSys(cs), ...
                        ControlSurf, 'Unif', false)';                    
                end
            end                 
            
            %Collect all the aerodynamic panel objects
            AeroPanels = [AeroPanel, ControlPanels];
            nAeroSeg = numel(AeroPanels);  
            
            %Assign a reference to the aerodynamic properties
            set(AeroPanels, 'AeroProp', AeroProp);
            
            %Preallocate
            StrucSet = arrayfun(@(~) awi.fe.StructuralSet    , 1 : nAeroSeg, 'Unif', false);   
            AeroSet  = arrayfun(@(~) awi.fe.AeroPanelSet     , 1 : nAeroSeg, 'Unif', false);
            Spline   = arrayfun(@(~) awi.fe.AeroelasticSpline, 1 : nAeroSeg, 'Unif', false); 
            StrucSet = horzcat(StrucSet{:});
            AeroSet  = horzcat(AeroSet{:});
            Spline   = horzcat(Spline{:});  
            
            %Set up coordinate matrices - Quicker outside the loop
            %   - Need absolute coordinate of beam coordinate and panel LE
            %     coordinates.
            aBeam = abs(beamX);
            eBeam = (beamX(ind, :) - beamX(ind, 1)) ./ obj.Span;
            aX1   = abs(horzcat(AeroPanels.X1));
            aX4   = abs(horzcat(AeroPanels.X4));
            
            %Must create new coordinate systems for the splines
            %   - These coordinate systems must have their y-axis along the
            %     beam axis.
            [cs, eta] = createCoordSysObjects(obj, 'GridVector', 'y');
            SplineCoordSys = makeFECoordSys(obj, cs(1 : end - 1));
            coordSysIndex  = 1 : numel(cs);
            
            %Loop through the aerodynamic panel sets and establish which
            %structural DOFs are required to construct the splines
            for iAP = 1 : numel(AeroPanels)
            
                %Index the coordinates within X1/X4 of this aero panel set
                idx = and(aBeam(ind, :) >= aX1(ind, iAP), aBeam(ind, :) <= aX4(ind, iAP));
               
                %Grab the nodes and add to the structural set
                beamNode = BeamNodes(idx);
                StrucSet(iAP).Nodes = horzcat(beamNode, NodesLE(idx), NodesTE(idx))';                
                
                %Select coordinate system which defines the beam axis
                %   - Each CAERO1 segment must be planar, therefore, each
                %     node in the 'beamNode' variable will have the same
                %     coordinate system defining the orientation of the
                %     beam. -> Can just reference the output coordinate
                %     system of the node as each node has its output
                %     coordinate system aligned with the beam axis.
                index = interp1(eta, coordSysIndex, eBeam(idx), 'previous');
                Spline(iAP).CoordSys = SplineCoordSys(index(1));
                
            end
            
            %Populate the aerodynamic cards
            %   - Each CAERO1 card has an AELIST with all the ID numbers
            %   - Each spline points to its corresponding SET1 entry
            %   - Each Splint points to its corresponding AELIST entry
            %     i.e. 1 x spline per CAERO1
            set(AeroSet, {'AeroPanels'}  , num2cell(AeroPanels)');
            set(Spline , {'AeroPanelSet'}, num2cell(AeroSet)'); 
            set(Spline , {'StrucSet'}    , num2cell(StrucSet)');            
            
            %Create the specific control surface cards
            if ~isempty(ControlPanels)
                nCS = numel(ControlPanels);
                %Make the panels
                ControlSurfs = arrayfun(@(~) awi.fe.AeroControlSurf, ...
                    1 : nCS, 'Unif', false);
                ControlSurfs = horzcat(ControlSurfs{:});
                %Assign the correct label
                set(ControlSurfs, {'LABEL'}, {ControlSurf.Label}');
                %Assign the aero panel sets
                set(ControlSurfs, {'AeroPanelSet'}, num2cell(AeroSet(end - (nCS-1) : end))');
                %Assign the coordinate systems
                set(ControlSurfs, {'CoordSys'}, hingeCoordSys');
                addFEData(FEModel, hingeCoordSys{:});
            end
        
            %Add objects to the FE model
            % - Add the aero control surfaces seperately to other panels
            addFEData(FEModel, AeroPanels, ControlSurfs, AeroProp, ...
                AeroSet, StrucSet, Spline, SplineCoordSys);
            
            return
        end
        %2D wingbox model
        function FEModel = generate2DFEModel(obj, FEModel, varargin)
            %generate1DBeamModel Generates the 'awi.fe.' data necessary to
            %describe a 2D panel element representation of a LiftingSurface
            %object based on the specified Spar/Stringer/Rib layout.
            
            type = obj.MeshType;
            assert(strcmp(type, 'uniform'), ['The ''MeshType'' must ', ...
                'be set to ''uniform''. Code update needed to '      , ...
                'handle generic unstructured mesh types.']);
            
            FEModel.log(sprintf('Beginning GFEM generation for object ''%s''...', obj.Name));
            
            %Debugging - Helps to see what the mesh looks like!
            bPlotMesh = true;
            
            [Spars, Ribs, Stringers, sw, sh] = i_parse(obj);
             
            function [Spars, Ribs, Stringers, sw, sh] = i_parse(obj)
                
                %Mesh size
                sw = obj.ShellWidth;
                ar = obj.ShellAR;
                if isempty(sw)
                    sw = 0.1;
                end
                if isempty(ar)
                    ar = 1;
                end
                sh = ar * sw;
                
                %Have the necessary components been defined?
                Spars     = find(obj.Children, 'Type', 'Spar');
                Stringers = find(obj.Children, 'Type', 'Stringer');
                Ribs      = find(obj.Children, 'Type', 'Rib');
                if isempty(Spars)
                    return
                end
                
                %Have the spars got enough information?
                etaThick = get(Spars, {'EtaThickness'});
                thickness = get(Spars, {'Thickness'});
                nData = [cellfun(@numel, etaThick), cellfun(@numel, thickness)];
                assert(all(all(nData > 0)), ['Unable ', ...
                    'to generate a GFEM as some of the Spars are '    , ...
                    'missing thickness or eta information.']);
                assert(all(range(nData) == 0), ['Unable to '   , ...
                    'generate a GFEM as some of the Spars do not have ' , ...
                    'the same number of points in the ''EtaThickness'' ', ...
                    'and ''Thickness'' properties. Update the '         , ...
                    'properties and try again.']);
                
            end
            
            if isempty(Spars)
                warning('Unable to generate GFEM model without spar objects.');
                return
            end
            if isempty(Ribs) || isempty(Stringers)
                FEModel.log(['No Ribs or Stringers defined - Using rib ', ...
                    'pitch and/or stringer pitch to define structural layout.']);
                %Use 'RibPitch' & 'StringPitch' to define wing box
                [Ribs, Stringers] = defineStructuralLayout(obj, Spars);
            end
               
            %% Define spar mesh using 'Ribs' geometry data
            
            FEModel.log('Meshing spars...');
            FEModel.HierarchyLevel = FEModel.HierarchyLevel + 1;
            
            [nBay, nEta]  = size(Ribs);
            nSp = nBay + 1;
            
            %Spar corner coordinates
            FEModel.log('Grabbing spar vertex positions');
            [xSp, ySp, zSp] = i_getSparCoords(Ribs, nBay, nEta, nSp);
            
            function [xSp, ySp, zSp] = i_getSparCoords(Ribs, nBay, nEta, nSp)
                %i_getSparCoords Returns the corner coordinates of the
                %rib-spar intersection points.
                %
                % Detailed Description:
                %   - 'Ribs' is a matrix of 'awi.model.CrossSection'
                %     objects. 
                %   - We know that the cross-sections will have a
                %     consistent number of points for the upper and lower
                %     points - this allows us to index into the coordinate
                %     data to find the coordinates at the LE & TE of the
                %     upper and lower surfaces of the rib.
                
                xSp = zeros(nSp, nEta, 2);
                ySp = zeros(nSp, nEta, 2);
                zSp = zeros(nSp, nEta, 2);
                for iBay = 1 : nBay
                    [xRib, yRib, zRib] = calculateGlobalCoords(Ribs(iBay, :));
                    nPoints = (size(xRib, 2) - 1) / 2;
                    %Upper surface
                    xSp(iBay    , :, 1) = xRib(:, 1);
                    ySp(iBay    , :, 1) = yRib(:, 1);
                    zSp(iBay    , :, 1) = zRib(:, 1);
                    xSp(iBay + 1, :, 1) = xRib(:, nPoints);
                    ySp(iBay + 1, :, 1) = yRib(:, nPoints);
                    zSp(iBay + 1, :, 1) = zRib(:, nPoints);
                    %Lower surface
                    xSp(iBay + 1, :, 2) = xRib(:, nPoints + 1);
                    ySp(iBay + 1, :, 2) = yRib(:, nPoints + 1);
                    zSp(iBay + 1, :, 2) = zRib(:, nPoints + 1);
                    xSp(iBay    , :, 2) = xRib(:, end - 1);
                    ySp(iBay    , :, 2) = yRib(:, end - 1);
                    zSp(iBay    , :, 2) = zRib(:, end - 1);
                end
                
            end
            
            %Vertical mesh seed
            FEModel.log('Populating spar vertical mesh seed');
            [meshCoords, nShellHeight] = i_getVerticalMeshSeed(xSp, ySp, zSp, nSp, nEta, sh, type);
            
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

            function [r, r_mod, r_norm] = i_getVector(x, y, z, dimDiff, dimCat)
                %i_getVector Calculates the vector ('r') between a set of
                %points. Also returns the length of the straight line
                %('r_mod') and the normalised direction vector ('r_norm').
                r      = cat(dimCat, diff(x, [], dimDiff), diff(y, [], dimDiff), diff(z, [], dimDiff));
                r_mod  = vecnorm(r, 2, dimCat);
                r_norm = r ./ r_mod;                
            end
            
            %Spanwise mesh (based on vertical mesh seed)
            FEModel.log('Calculating spar spanwise nodes based on vertical mesh seed');
            sparMesh = i_getSparMesh(xSp, ySp, zSp, meshCoords, nSp, nEta, nShellHeight, sw, type);
           
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

            %Plot to check?
            if bPlotMesh
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

            FEModel.HierarchyLevel = FEModel.HierarchyLevel - 1;
            
            %% Define the rib mesh 
            
            FEModel.log('Meshing ribs...');
            FEModel.HierarchyLevel = FEModel.HierarchyLevel + 1;
            
            assert(strcmp(type, 'uniform'), ['The ''MeshType'' must ', ...
                'be set to ''uniform''. Code update needed to '      , ...
                'handle generic unstructured mesh types.']);
            
            %Get spar-rib intersection points and calculate no. shells
            FEModel.log('Finding rib/spar intersection points');
            [xRibSpar, yRibSpar, zRibSpar, nShellRib] = i_getRibCornerCoords(sparMesh, xSp, sw);
            
            function [xRibSpar, yRibSpar, zRibSpar, nShellRib] = i_getRibCornerCoords(sparMesh, xSp, sw)
                
                %Where do the ribs intersect the spars?
                indRib{1} = find(ismember(sparMesh{1}(1, :, 1), xSp(1, :, 1)));
                indRib{2} = find(ismember(sparMesh{2}(1, :, 1), xSp(2, :, 1)));
                indRib{3} = find(ismember(sparMesh{3}(1, :, 1), xSp(3, :, 1)));
                assert(range(cellfun(@numel, indRib)) == 0, ['Expected there ', ...
                    'to be the same number of rib intersections for each spar']);
                indRib = vertcat(indRib{:});
                assert(nnz(range(indRib)) == 0, ['Expected the ribs to ' , ...
                    'intersect each spar at the same spanwise position. ', ...
                    'Make sure the spar mesh is uniform.']);
                indRib = indRib(1, :);
                
                %Extract spar corner points
                sparCorners = cellfun(@(x) x(:, indRib, :), sparMesh, 'Unif', false);
                sparCorners = cat(4, sparCorners{:});
                xRibSpar = squeeze(sparCorners(:, :, 1, :));
                yRibSpar = squeeze(sparCorners(:, :, 2, :));
                zRibSpar = squeeze(sparCorners(:, :, 3, :));
                
                %Calculate number of panels per rib
                [~, r_mod, ~] = i_getVector(xRibSpar, yRibSpar, zRibSpar, 3, 4);
                nShellRib = squeeze(max(max(ceil(r_mod ./ sw), [], 2), [], 1));
                
            end
                  
            %Define coordinates of nodes on each rib face
            FEModel.log('Calculating rib mesh');
            ribMesh = i_calculateRibMesh(Ribs, xRibSpar, yRibSpar, zRibSpar, nShellRib, nShellHeight, nBay);
            
            function ribMesh = i_calculateRibMesh(Ribs, xRibSpar, yRibSpar, zRibSpar, nShellRib, nShellHeight, nBay)
                
                nRib = size(Ribs, 2);
                
                %For each bay...
                %   - Interpolate the spar profile at desired panel positions
                %   - Generate new nodes along the face of each rib
                %   - Do NOT generate duplicate nodes at the spar/rib locs
                ribMesh = arrayfun(@(i) zeros((nShellRib(i) - 1) * (nShellHeight(i) + 1), nRib, 3), 1 : nBay, 'Unif', false);
                for iBay = 1 : nBay
                    
                    %Normalise the profile
                    [xNorm, zNorm, xzData] = normaliseXZData(Ribs(iBay, :));
                    
                    %Split into upper and lower surfaces
                    nPoints = (size(xNorm, 2) - 1) / 2;
                    xU = xNorm(:, 1 : nPoints);
                    zU = zNorm(:, 1 : nPoints);
                    xL = fliplr(xNorm(:, nPoints + 1 : end - 1));
                    zL = fliplr(zNorm(:, nPoints + 1 : end - 1));
                    
                    %Interpolate the upper and lower surfaces
                    xRibU = repmat(linspace(0, 1, nShellRib(iBay) + 1), [nRib, 1]);
                    xRibL = xRibU;
                    zRibU = zeros(nRib, nShellRib(iBay) + 1);
                    zRibL = zRibU;
                    for iSpan = 1 : size(xU, 1)
                        zRibU(iSpan, :) = interp1(xU(iSpan, :), zU(iSpan, :), xRibU(iSpan, :));
                        zRibL(iSpan, :) = interp1(xL(iSpan, :), zL(iSpan, :), xRibL(iSpan, :));
                    end
                    
                    %Scale/translato the coords back to correct chord/position
                    [xRibU, zRibU] = scaleProfileData(xRibU, zRibU, xzData);
                    [xRibL, zRibL] = scaleProfileData(xRibL, zRibL, xzData);
                    
                    %Assign back to the cross-sections so we can use methods
                    set(Ribs(iBay, :), {'X'}, num2cell([xRibU , fliplr(xRibL)], 2));
                    set(Ribs(iBay, :), {'Z'}, num2cell([zRibU , fliplr(zRibL)], 2));
                    [xRib, yRib, zRib] = calculateGlobalCoords(Ribs(iBay, :));
                    
                    %Get directon vector between each node on the upper and
                    %lower surface
                    nPoints = size(xRib, 2) / 2;
                    xRib = cat(3, xRib(:, 1 : nPoints), fliplr(xRib(:, nPoints + 1 : end)));
                    yRib = cat(3, yRib(:, 1 : nPoints), fliplr(yRib(:, nPoints + 1 : end)));
                    zRib = cat(3, zRib(:, 1 : nPoints), fliplr(zRib(:, nPoints + 1 : end)));
                    [~, r_mod, r_norm] = i_getVector(xRib, yRib, zRib, 3, 3);
                    
                    %Define new nodes for each rib using the direction
                    %vector (vectorised)
                    %ub = cumsum(repmat(nShellHeight(iBay) + 1, [1, nShellRib(iBay)]));
                    %lb = [1, ub(1 : end - 1) + 1];
                    %    - First set of nodes will be the same as the spar
                    %ribMesh{iBay}(lb(1)   : ub(1)  , :, :) = cat(3, xRibSpar(:, :, iBay)    , yRibSpar(:, :, iBay)    , zRibSpar(:, :, iBay));
                    %ribMesh{iBay}(lb(end) : ub(end), :, :) = cat(3, xRibSpar(:, :, iBay + 1), yRibSpar(:, :, iBay + 1), zRibSpar(:, :, iBay + 1));
                    %    - Equally spaced nodes along rib
                    ds    = permute(linspace(0, 1, nShellHeight(iBay) + 1), [1, 3, 2]);
                    dR    = r_mod(:, 2 : end - 1) .* ds;                          
                    %   - Calculate node positions
                    ribMesh{iBay}(:, :, 1) = i_calculateRibFaceNodes(xRib(:, 2 : end - 1, 1), r_norm(:, 2 : end - 1, 1), dR);
                    ribMesh{iBay}(:, :, 2) = i_calculateRibFaceNodes(yRib(:, 2 : end - 1, 1), r_norm(:, 2 : end - 1, 2), dR);
                    ribMesh{iBay}(:, :, 3) = i_calculateRibFaceNodes(zRib(:, 2 : end - 1, 1), r_norm(:, 2 : end - 1, 3), dR);
                    
                end
                
                %Combine the coordinates from each bay to form one single
                %rib!
                %endSpar = cat(3, xRibSpar(:, :, end), yRibSpar(:, :, end), zRibSpar(:, :, end));
                %ribMesh = vertcat(ribMesh{:}, endSpar);
            end
            
            function coords = i_calculateRibFaceNodes(r0, r_norm, dR)
                %i_calculateRibFaceNodes Calculates the x, y or z position
                %of the rib face nodes using the direction vector of the
                %rib face and the desired rib node spacing.
                
                %r = r0 + ds * r_norm
                coords = num2cell(permute( ...
                    r0 + r_norm .* dR, [3, 2, 1]), [1, 2]);
                coords = cellfun(@(x) x(:), coords, 'Unif', false);
                coords = horzcat(coords{:});
                
            end
            
            if bPlotMesh
                ribData = vertcat(ribMesh{:});
                xRib_ = i_flatten(ribData(:, :, 1));
                yRib_ = i_flatten(ribData(:, :, 2));
                zRib_ = i_flatten(ribData(:, :, 3));
                plot3(xRib_, yRib_, zRib_, ...
                    'Marker', 'o', 'MarkerEdgeColor', 'k', ...
                    'MarkerFaceColor', 'r', 'LineStyle', 'none', ...
                    'DisplayName', 'Rib Face Nodes');
            end
            
            FEModel.HierarchyLevel = FEModel.HierarchyLevel - 1;
            
            %% Define the skin mesh
            
            %Make new 'Rib' cross sections at every location where the spar
            %has a node --> Need to interpolate the profiles.
            
            %Mesh in a similar fashion to the rib object
            
            %% Generate 'awi.fe' objects
            % Things to do:
            %   1. Define panel properties
            %   2. Define panel materials
            %   3. Generate RBE spider at the inboard rib/apply SPCs
            %   4. Define skin mesh             
            
            FEModel.log('Instantiating FE objects...');
            FEModel.HierarchyLevel = FEModel.HierarchyLevel + 1;
            
            %Spars
            FEModel.log('Spars...');
            [SparNodes, SparPanels, SparProps] = i_makeSparFE(sparMesh, Spars);
            
            function [SparNodes, SparPanels, SparProps] = i_makeSparFE(sparMesh, Spars)
                
                SparProps = awi.fe.PanelPropCollection;
                
                %How many nodes & panels?
                nNodePerSpar  = cellfun(@(x) numel(x(:, :, 1)), sparMesh);
                nPanelPerSpar = cellfun(@(x) (size(x(:, :, 1), 1) - 1) * (size(x(:, :, 1), 2) - 1), sparMesh);
                nSparNodes    = sum(nNodePerSpar);
                nSparPanels   = sum(nPanelPerSpar);
                
                %Indexing
                bound = [0, cumsum(nNodePerSpar(1 : end - 1))];
                ubSp  = cumsum(nNodePerSpar);
                lbSp  = [1, ubSp(1 : end - 1) + 1];
                ubPa  = cumsum(nPanelPerSpar);
                lbPa  = [1, ubPa(1 : end - 1) + 1];
                
                %Make the objects (1 collection)
                SparNodes   = awi.fe.NodeCollection;
                SparNodes.X = zeros(3, nSparNodes);
                SparPanels  = awi.fe.PanelCollection;
                SparPanels.NodeCollection = SparNodes;
                SparPanels.GridID         = zeros(4, nSparPanels);
                
                %For each spar...
                for iiSp = 1 : numel(sparMesh)
                    
                    %Format coords for each spar
                    [X, index] = i_stackCoords(sparMesh{iiSp});
                    SparNodes.X(:, lbSp(iiSp) : ubSp(iiSp)) = X;
                    
                    %Define the shells
                    %   - Index into the node collection
                    g1 = index(1 : end - 1, 1 : end - 1);
                    g2 = index(1 : end - 1, 2 : end);
                    g3 = index(2 : end    , 2 : end);
                    g4 = index(2 : end    , 1 : end - 1);
                    nodeIndex = [g1(:), g2(:), g3(:), g4(:)]' + bound(iiSp);
                    %   - Assign data
                    SparPanels.NodeIndex(:, lbPa(iiSp) : ubPa(iiSp)) = nodeIndex;
                    
                    %Add a subgroup
                    str = sprintf('Spar %i', iiSp);
                    SparPanels.SubGroups = [SparPanels.SubGroups, {str, nodeIndex}];
                    
                    %Define shell properties
                    
                end
                
            end
                          
            function [coords, index] = i_stackCoords(mesh)
                [nr, nc, ~] = size(mesh);
                x_ = num2cell(mesh(:, :, 1), 1);
                y_ = num2cell(mesh(:, :, 2), 1);
                z_ = num2cell(mesh(:, :, 3), 1);
                coords = [vertcat(x_{:}), vertcat(y_{:}), vertcat(z_{:})]';
                index = reshape(1 : (nr * nc), [nr, nc]);
            end
            
            %Ribs
            FEModel.log('Ribs...');
            [RibAndSparNodes, RibPanels, SparPanels] = i_makeRibFE( ...
                SparNodes, SparPanels, xRibSpar, yRibSpar, zRibSpar, ribMesh, nShellHeight, nShellRib, nBay);
             
            %[RibNodes, RibPanels] = i_makeRibFE(ribMesh, nShellRib, nShellHeight, nBay);
            
            function [RibAndSparNodes, RibPanels, SparPanels] = i_makeRibFE(SparNodes, SparPanels, xRibSpar, yRibSpar, zRibSpar, ribMesh, nShellHeight, nShellRib, nBay)
                
                sparCoords = SparNodes.X';
                
                nSparNodes = size(SparNodes.X, 2);
                nRib       = size(xRibSpar, 2);
                
                %Make a new collection for the combined spar & rib elements
                RibAndSparNodes          = awi.fe.NodeCollection;
                RibPanels                = awi.fe.PanelCollection;
                RibPanels.NodeCollection = RibAndSparNodes;
                
                %Reassign the collection for the Spar panels
                SparPanels.NodeCollection = RibAndSparNodes;
                
                %How many Nodes/Panels per rib section?
                nNodePerRib    = cellfun(@(x) size(x, 1), ribMesh);
                nNodePerRibSet = cellfun(@(x) numel(x(:, :, 1)), ribMesh);
                nPanelPerRib   = arrayfun(@(i) nShellRib(i) * nShellHeight(i), 1 : numel(ribMesh));
                
                %Combine the rib coordinates and add to the collection
                ribCoords = arrayfun(@(i)  [ ...
                    reshape(ribMesh{i}(:, :, 1), [1, nNodePerRibSet(i)]) ; ...
                    reshape(ribMesh{i}(:, :, 2), [1, nNodePerRibSet(i)]) ; ...
                    reshape(ribMesh{i}(:, :, 3), [1, nNodePerRibSet(i)])], 1 : nBay, 'Unif', false);
                RibAndSparNodes.X = [SparNodes.X, horzcat(ribCoords{:})];
                
                %Indexing
                bayBounds = [0, cumsum(nNodePerRibSet(1 : end - 1))] + nSparNodes;                
                
                %Loop through each bay & rib and assign index number for Rib/Spar node
                %collection
                %   - Need to use a loop as we need to find the indices of the rib/spar
                %     connection points for each rib. Not sure how to vectorise...
                nodeIndex = arrayfun(@(n) repmat({zeros(4, n)}, [1, nRib]), nPanelPerRib, 'Unif', false);
                nodeIndex = vertcat(nodeIndex{:});
                for iB = 1 : nBay
                    
                    %Indexing
                    ribBounds = [0, cumsum(repmat(nNodePerRib(iB), [1, nRib - 1]))];
                    ind0      = reshape((1 : nNodePerRib(iB)), [nShellHeight(iB) + 1, nShellRib(iB) - 1]);
                    index0    = arrayfun(@(n) ind0 + n, ribBounds, 'Unif', false);
                    
                    for iR = 1 : nRib
                        
                        %Coordinates of LE/TE of rib
                        ribLE = [xRibSpar(:, iR, iB)    , yRibSpar(:, iR, iB)    , zRibSpar(:, iR, iB)];
                        ribTE = [xRibSpar(:, iR, iB + 1), yRibSpar(:, iR, iB + 1), zRibSpar(:, iR, iB + 1)];
                        
                        %Find the LE/TE in the spar coordinates
                        indLE = find(ismember(sparCoords, ribLE, 'rows'));
                        indTE = find(ismember(sparCoords, ribTE, 'rows'));
                        
                        assert(and(any(indLE), any(indTE)), sprintf([ ...
                            'Unable to find the LE/TE connection points for rib no. %i ', ...
                            'in bay no. %i. Check mesh.'], iR, iB));
                        
                        %Add the index numbers of nodes on the rib LE/TE
                        index = [indLE, index0{iR} + bayBounds(iB), indTE];
                        
                        %Assign panel vertices
                        g1  = index(1 : end - 1, 1 : end - 1);
                        g2  = index(1 : end - 1, 2 : end);
                        g3  = index(2 : end    , 2 : end);
                        g4  = index(2 : end    , 1 : end - 1);
                        nodeIndex{iB, iR} = [g1(:), g2(:), g3(:), g4(:)]';
                        
                    end
                end
                
                %Define SubGroups?
                RibPanels.SubGroups = {'Inboard Rib', horzcat(nodeIndex{:, 1})};
                
                %Assign to PanelCollection object
                allRibNodeIndex = arrayfun(@(i) horzcat(nodeIndex{i, :}), 1 : nBay, 'Unif', false);
                RibPanels.NodeIndex = horzcat(allRibNodeIndex{:});
                
            end
            
            function [RibNodes, RibPanels] = i_makeRibFE_OLD(ribMesh, nShellRib, nShellHeight, nBay)
                
                %This old code will is applicable for the case where the
                %'ribMesh' has duplicate nodes where the rib/spar connects.
                %This could be useful if we ever want to model fasteners
                %(possibly using CBUSH elements). Perhaps ask Lucien about
                %this.
                
                %nRib = size(ribMesh{1}, 2);
                nRib = size(ribMesh, 2);
                
                assert(and(range(nShellRib) == 0, range(nShellHeight) == 0), ...
                    ['Expected the rib mesh to be uniform. Update code ', ...
                    'for unstructured mesh.']);
                
                %How many nodes & panels?
                nNodesPerRib = size(ribMesh, 1);
                nRibNodes    = numel(ribMesh(:, :, 1));
                nChordPanel  = ((nNodesPerRib / (nShellHeight(1) + 1)) - 1);
                nPanelPerRib = nChordPanel * nShellHeight(1);
                nRibPanels   = nPanelPerRib * nRib;
                %nRibPanels = (size(ribMesh(:, :, 1), 1) - 1) * (size(ribMesh(:, :, 1), 2) - 1);
                %nNodePerRib     = cellfun(@(x) numel(x(:, :, 1)), ribMesh);
                %nPanelPerRibSet = arrayfun(@(i) nShellRib(i) * nShellHeight(i) * size(ribMesh{i}, 2), 1 : numel(ribMesh));
                %nRibNodes   = sum(nNodePerRib);
                %nRibPanels  = sum(nPanelPerRibSet);
                
                %Indexing
                bound = [0, cumsum(repmat(nNodesPerRib, [1, nRib - 1]))];
%                 bound = [0, cumsum(nNodePerRib(1 : end - 1))];
%                 ubSp  = cumsum(nNodePerRib);
%                 lbSp  = [1, ubSp(1 : end - 1) + 1];
%                 ubPa  = cumsum(nPanelPerRibSet);
%                 lbPa  = [1, ubPa(1 : end - 1) + 1];
                
                %Make the objects (1 collection)
                RibNodes   = awi.fe.NodeCollection;
                RibPanels  = awi.fe.PanelCollection;
                RibPanels.NodeCollection = RibNodes;
                RibPanels.GridID         = zeros(4, nRibPanels);
                
                %Stack the coordinates
                RibNodes.X = [ ...
                    reshape(ribMesh(:, :, 1), [1, nRibNodes]) ; ...
                    reshape(ribMesh(:, :, 2), [1, nRibNodes]) ; ...
                    reshape(ribMesh(:, :, 3), [1, nRibNodes])];                
                                
                %Define indexing for the rib panels
                ind = reshape((1 : nNodesPerRib), [nShellHeight(1) + 1, nChordPanel + 1]);
                g1  = ind(1 : end - 1, 1 : end - 1);
                g2  = ind(1 : end - 1, 2 : end);
                g3  = ind(2 : end    , 2 : end);
                g4  = ind(2 : end    , 1 : end - 1);
                nodeIndex = [g1(:), g2(:), g3(:), g4(:)]';
                %   - Increase index number by nNodesPerRib for each rib
                nodeIndex = arrayfun(@(n) nodeIndex + n, bound, 'Unif', false);
                RibPanels.NodeIndex = horzcat(nodeIndex{:});
                
%                 %For each chordwise bay                
%                 for iBay = 1 : nBay
%                     nH = nShellHeight(iBay) + 1;
%                     nR = nShellRib(iBay) + 1;
%                     
%                     %Format coords for each spar
%                     [X, index] = i_stackCoords(ribMesh{iBay});
%                     RibNodes.X(:, lbSp(iBay) : ubSp(iBay)) = X;
%                     
%                     %Index each rib
%                     nodeIndex = num2cell(zeros(4, nShellRib(iBay) * nShellHeight(iBay), nRib), [1, 2]);
%                     for i = 1 : nRib
%                         
%                         ind = reshape(index(:, i), [nH, nR]);
%                                                 
%                         %Define the shells
%                         %   - Index into the node collection
%                         g1 = ind(1 : end - 1, 1 : end - 1);
%                         g2 = ind(1 : end - 1, 2 : end);
%                         g3 = ind(2 : end    , 2 : end);
%                         g4 = ind(2 : end    , 1 : end - 1);
%                         
%                         nodeIndex{i} = [g1(:), g2(:), g3(:), g4(:)]';
%                         
%                     end                    
%                     nodeIndex = horzcat(nodeIndex{:}) + bound(iBay);
%                     
%                     %Assign data
%                     RibPanels.NodeIndex(:, lbPa(iBay) : ubPa(iBay)) = nodeIndex;    
%                     
%                 end
                
                
            end            
                        
            %Plot to check
            if bPlotMesh
                hF  = figure;
                hAx = axes('Parent', hF);
                drawElement(RibAndSparNodes, hAx);
                drawElement(SparPanels, hAx);
                drawElement(RibPanels , hAx);
            end          
            
            FEModel.HierarchyLevel = FEModel.HierarchyLevel - 1;
            
        end
        %Aero panel methods
        function AeroPanel = defineAeroPanelsNoControlSurf(obj, PanelCoord)
            %i_defineAeroPanelsNoControlSurf Defines the aerodynamic
            %panels for the case where the lifting surface does not
            %have any control surfaces.
            
            %Better if the panel coordinates are passed in
            if nargin < 2
               PanelCoord = obj.PanelCoords; 
            end
            
            %Get coords of LE & TE (remember, this is in the global frame!)
            %   - TODO : What is this is empty??            
            LE_ = [PanelCoord.LE.X ; PanelCoord.LE.Y ; PanelCoord.LE.Z];
            TE_ = [PanelCoord.TE.X ; PanelCoord.TE.Y ; PanelCoord.TE.Z];
            
            %How many aerodynamic sections are there?
            nAeroSeg = size(LE_, 2) - 1;
            
            %Determine which set of coordinate to use when indexing the
            %beam coordinates
            switch obj.SpanVector
                case 'X'
                    error('Unable to proceed for SpanVector = ''X''. Update code.');
                case 'Y'
                    index = 2;
                case 'Z'
                    index = 3;
            end
            
            %Define initial object
            AeroPanel = arrayfun(@(~) awi.fe.AeroPanel, 1 : nAeroSeg, 'Unif', false);
            AeroPanel = horzcat(AeroPanel{:});
            
            %Determine segment length in spanwise and chordwise direction
            s = abs(diff(LE_(index, :), 1));
            c = abs(LE_(1, :) - TE_(1, :));
            chrd = [c(1 : end - 1) ; c(2 : end)];
            c = max(chrd, [], 1);
            [nSpan, nChord] = obj.getAeroPanelSize(obj, s, c);
            
            %Determine the twist distribution along each panel segment
            eta_ = unique([obj.Eta_,obj.Twist_eta]);
            twist_ = interp1(obj.Twist_eta,obj.Twist,eta_);
            for i = 1:nAeroSeg
                I = eta_>=obj.Eta_(i) & eta_<=obj.Eta_(i+1);
                tmp_twist = twist_(I);
                tmp_eta = eta_(I);
                tmp_eta = (tmp_eta - min(tmp_eta))./(max(tmp_eta)-min(tmp_eta));
                AeroPanel(i).AoA = tmp_twist;
                AeroPanel(i).AoA_eta = tmp_eta;
            end
            
            %Populate the 'AeroPanels' structure
            set(AeroPanel, {'NSPAN'} , num2cell(nSpan)');
            set(AeroPanel, {'NCHORD'}, num2cell(nChord)');
            set(AeroPanel, {'X1'}    , num2cell(LE_(:, 1 : end - 1), 1)');
            set(AeroPanel, {'X4'}    , num2cell(LE_(:, 2 : end), 1)');
            set(AeroPanel, {'CHORD'} , num2cell(chrd, 1)');
            
        end        
        function AeroPanel = defineAeroPanelsWithControlSurf(obj, ControlSurf, PanelCoord)
            %i_defineAeroPanelsWithControlSurf Defines the aerodynamic
            %panels for the case where we have control surfaces.
            
            %Better if the panel coordinates are passed in
            if nargin < 3
               PanelCoord = obj.PanelCoords; 
            end
            
            %Parse the control surfaces
            Dev = checkControlSurfaceOverlap(obj, ControlSurf);
            
            %Sort all devices from root to tip
            devEta = vertcat(Dev.Eta);
            [~, ind] = sort(devEta(:, 1));
            Dev = Dev(ind);
            devEta = vertcat(Dev.Eta);
            
            %There will need to be a new panel at every unique eta position
            eta  = unique([devEta(:) ; [0 ; 1] ; obj.Chord_eta(:)]);
            %Grab chord value and global position of LE at each eta
            chrd = interp1(obj.Chord_eta, obj.Chord, eta);
            LEx  = interp1(obj.Eta_, PanelCoord.LE.X, eta);
            
            %Determine which control surfaces are present at the spanwise
            %locations
            e1 = eta(1 : end - 1)';
            eN = eta(2 : end)';
            idx = or( ...
                and(e1 >= devEta(:, 1), e1 <= devEta(:, 2)), ...
                and(eN >= devEta(:, 1), eN <= devEta(:, 2)));
            dev = arrayfun(@(i) Dev(idx(:, i)), 1 :  numel(eta) - 1, 'Unif', false);
            
            %Determine the absolute x-coordinate of each corner of each
            %panel segment and the chord length of each side of the panel
            xLE   = cell(1, numel(dev));
            chord = cell(1, numel(dev));
            for ii = 1 : numel(dev) %Loop through eta positions
                if isempty(dev{ii})
                    %Assume panels run for the full-chord length at the
                    %inboard and outboard positions
                    xle = [0, 0];
                    xte = [1, 1];
                else
                    %Grab the start and end positions of the segment
                    e = [eta(ii), eta(ii + 1)];
                    %Of the control surfaces that neighbour this segment,
                    %which ones sit either-side?
                    devEta = vertcat(dev{ii}.Eta);
                    isNeighbour = or( ...
                        and(devEta(:, 1) == e(2), devEta(:, 2) > e(2)), ...
                        and(devEta(:, 2) == e(1), devEta(:, 1) < e(1)));
                    %Down-select control surfaces to remove neighbours
                    cs = dev{ii}(~isNeighbour);
                    %Determine the normalised offset to the LE and TE for
                    %each control surface at the current eta position.
                    xle_ = arrayfun(@(d) interp1(d.Eta, d.xLE, e, 'linear', 'extrap'), cs, 'Unif', false);
                    xte_ = arrayfun(@(d) interp1(d.Eta, d.xTE, e, 'linear', 'extrap'), cs, 'Unif', false);
                    xle_ = vertcat(xle_{:});
                    xte_ = vertcat(xte_{:});
                    %Find unique quartets - Can sometimes get duplicates
                    [~, ind] = uniquetol([xle_, xte_], 1e-4, 'ByRows', true);
                    xle_ = xle_(ind, :);
                    xte_ = xte_(ind, :);
                    %Store the 4 vertices of each control surface
                    cs_le_te = permute(cat(3, xle_, xte_), [3, 2, 1]);
                    %Add [0,1] to the list of positions - Means we never
                    %miss a part of the planform
                    xle = [[0,0] ; xte_];
                    xte = [xle_ ; [1, 1]];
                    %Sort the positions from LE -> TE
                    [~, ind] = sort(xle(:, 1), 'ascend');
                    xle = xle(ind, :);
                    xte = xte(ind, :);                    
                    %Create the LE/TE pairs for each segment
                    corners = [xle ; xte];
                    corners = uniquetol(corners, 1e-4, 'ByRows', true);                    
                    %Return to le/te format
                    le = corners(1 : end - 1, :);
                    te = corners(2 : end    , :);
                    %Store the vertices of each candidate panel segment
                    le_te = permute(cat(3, le, te), [3, 2, 1]);                    
                    %Compare against control surfaces and determine which LE/TE
                    %pairs belong to control surfaces and which are just
                    %lifting surfaces.
                    % - This is essentially two nested loops. We loop through
                    %   each segment and check if any of the control surfaces
                    %   have the same corner coordinates. If they do then we
                    %   remove them.
                    idx = arrayfun(@(iS) ...
                        any(arrayfun(@(iC) isequal( ...
                        le_te(:, :, iS), cs_le_te(:, :, iC)), ...
                        1 : size(cs_le_te, 3))), 1 : size(le_te, 3));
                    le_te_ = le_te(:, :, ~idx);
                    %Return to LE/TE format
                    xle = permute(le_te_(1, :, :), [3, 2, 1]);
                    xte = permute(le_te_(2, :, :), [3, 2, 1]);                    
                end
                %Calculate the x-coordinates of each panel segment at the
                %inboard and outboard positions and use this to calculate
                %the chord length
                ch   = repmat(chrd([ii, ii + 1])', [size(xle, 1), 1]);
                xOff = repmat(LEx([ii, ii+ 1])', [size(xle, 1), 1]);
                xLE{ii} = (xle .* ch) + xOff;
                xTE_ = (xte .* ch) + xOff;
                chord{ii} = diff(cat(3, xLE{ii}, xTE_), [], 3);
                %Check if any segments have been defined where the start
                %and end chord is equal to zero - Remove these rows
                idx = ismember(chord{ii}, [0, 0], 'rows');
                xLE{ii}   = xLE{ii}(~idx, :);
                chord{ii} = abs(chord{ii}(~idx, :));
            end
            
            %(y,z) coordinates are simple
            pos = obj.AbsPosition;
            switch obj.SpanVector
                case 'X'
                    error('Update code for the case where ''SpanVector'' is ''X''.');
                case 'Y'
                    yLE = eta .* obj.Span + pos(2);
                    zLE = awi.model.Component.interp1( ...
                        obj.YData ./ obj.Span, obj.ZData,  eta) + pos(3);
                case 'Z'
                    yLE = awi.model.Component.interp1( ...
                        obj.ZData ./ obj.Span, obj.YData,  eta) + pos(2);
                    zLE = eta .* obj.Span + pos(3);
            end
            
            %Assign the coordinates for each panel
            y = [yLE(1 : end - 1) , yLE(2 : end)];
            z = [zLE(1 : end - 1) , zLE(2 : end)];
            yLE = arrayfun(@(i) repmat(y(i, :), [size(xLE{i}, 1), 1]), ...
                1 : numel(xLE), 'Unif', false);
            zLE = arrayfun(@(i) repmat(z(i, :), [size(xLE{i}, 1), 1]), ...
                1 : numel(xLE), 'Unif', false);
            
            %Stack the coordinates
            xLE  = vertcat(xLE{:});
            yLE  = vertcat(yLE{:});
            zLE  = vertcat(zLE{:});
            chord = vertcat(chord{:});
            X1 = [xLE(:, 1), yLE(:, 1), zLE(:, 1)]';
            X4 = [xLE(:, 2), yLE(:, 2), zLE(:, 2)]';
            
            %Determine segment length in chordwise direction
            switch obj.SpanVector
                case 'X'
                    error('Update code for the case where ''SpanVector'' is ''X''.');
                case 'Y'
                    s = abs(diff(yLE, [], 2))';
                case 'Z'
                    s = abs(diff(zLE, [], 2))';
            end
            c = max(abs(chord), [], 2)';
            
            %Calculate aero panel size
            [nSpan, nChord] = obj.getAeroPanelSize(obj, s, c);
            
            %Generate the aerodynamic panels
            AeroPanel = arrayfun(@(~) awi.fe.AeroPanel, 1 : numel(nSpan), 'Unif', false);
            AeroPanel = horzcat(AeroPanel{:});
            
            %Determine the twist distribution along each panel segment
            eta_ = unique([eta(:);obj.Twist_eta(:)])';
            twist_ = interp1(obj.Twist_eta,obj.Twist,eta_);
            for i = 1:numel(nSpan)
                I = eta_>=eta(i) & eta_<=eta(i+1);
                tmp_twist = twist_(I);
                tmp_eta = eta_(I);
                tmp_eta = (tmp_eta - min(tmp_eta))./(max(tmp_eta)-min(tmp_eta));
                AeroPanel(i).AoA = tmp_twist;
                AeroPanel(i).AoA_eta = tmp_eta;
            end
            
            %Determine the AoA distribution along each panel segment
            %
            set(AeroPanel, {'NSPAN'} , num2cell(nSpan)');
            set(AeroPanel, {'NCHORD'}, num2cell(nChord)');
            set(AeroPanel, {'X1'}    , num2cell(X1, 1)');
            set(AeroPanel, {'X4'}    , num2cell(X4, 1)');
            set(AeroPanel, {'CHORD'} , num2cell(chord', 1)');
            
        end        
        function AeroPanel = defineControlSurfacePanels(obj, ControlSurf)
            %defineControlSurfacePanels Generates the 'awi.fe.AeroPanel'
            %objects for the control surfaces belonging to the lifting
            %surface.
            
            if nargin < 2
                ControlSurf = grabControlSurface(obj, ...
                    'awi.model.ControlSurface');
            end
            
            %Parse the control surfaces
            ControlSurf = checkControlSurfaceOverlap(obj, ControlSurf);
            
            %Grab the vertex coordinates in the global coordinate system
            %   - Let the 'awi.model.ControlSurface' object do all the
            %   heavy lifting
            %   - TODO: Vectorise this
            vertexCoords = [ControlSurf.Coords];
            
            %Check if any of the control surfaces have more than 2
            %coordinates in each direction
            idx = arrayfun(@(s) numel(s.LE.X) > 2, vertexCoords);
            
            %Those control surfaces with more coordinate points will need
            %to be split into multiple aero panels and then linked together
            linkSurf   = vertexCoords(idx);
            noLinkSurf = vertexCoords(~idx);
                        
            %Deal with the unlinked control surfaces first
%             AeroPanel = i_generatePanels(obj, noLinkSurf);

            if ~isempty(noLinkSurf) 
                 AeroPanel = i_generatePanels(obj, noLinkSurf);   
            end
            
            if isempty(linkSurf) %Escape route
                return
            end
            
            %Split each 'linkSurf' coord set into a standard coord set
            n = numel(linkSurf);
            Panel = cell(1, n);
            for ii = 1 : n
                nSeg = numel(linkSurf(ii).LE.X) - 1;
                le = arrayfun(@(i) struct( ...
                    'X', linkSurf(ii).LE.X([i, i + 1]) , ...
                    'Y', linkSurf(ii).LE.Y([i, i + 1]) , ...
                    'Z', linkSurf(ii).LE.Z([i, i + 1])), 1 : nSeg, 'Unif', false); 
                te = arrayfun(@(i) struct( ...
                    'X', linkSurf(ii).TE.X([i, i + 1]) , ...
                    'Y', linkSurf(ii).TE.Y([i, i + 1]) , ...
                    'Z', linkSurf(ii).TE.Z([i, i + 1])), 1 : nSeg, 'Unif', false);
                Panel{ii} = cell2struct([le ; te], {'LE', 'TE'});
            end
            
            %Make the panels
            AeroPanel_ = i_generatePanels(obj, vertcat(Panel{:})');
            
            %Return all panels 
            if ~isempty(noLinkSurf)
                AeroPanel = [AeroPanel, AeroPanel_];
            else
                AeroPanel = AeroPanel_;
            end
            
              function AeroPanel = i_generatePanels(obj, cs)
                %i_generatePanels Makes the 'awi.fe.AeroPanel' objects
                %using the LE & TE coordinate stored in the structure-array
                %'cs'.
                
                %Collect the LE and TE coordinates
                leCoords = [cs.LE];
                teCoords = [cs.TE];
                
                %Calculate coordinates of point (1) & point (4) and the chord
                %length at these points
                chrd = abs(vertcat(leCoords.X) - vertcat(teCoords.X));
                X1   = arrayfun(@(s) [s.LE.X(1) ; s.LE.Y(1) ; s.LE.Z(1)], cs, 'Unif', false);
                X4   = arrayfun(@(s) [s.LE.X(2) ; s.LE.Y(2) ; s.LE.Z(2)], cs, 'Unif', false);
                
                %Calculate the spanwise distance for each panel
                switch obj.SpanVector
                    case 'X'
                        error('Update code for ''SpanVector'' = X');
                    case 'Y'
                        s = abs(diff(vertcat(leCoords.Y), [], 2));
                    case 'Z'
                        s = abs(diff(vertcat(leCoords.Z), [], 2));
                end
                c = max(chrd, [], 2);
                
                %Calculate aero panel size
                [nSpan, nChord] = obj.getAeroPanelSize(obj.ControlSurfaces, s, c);
                
                %Generate the aerodynamic panels
                AeroPanel = arrayfun(@(~) awi.fe.AeroPanel, 1 : numel(cs), 'Unif', false);
                AeroPanel = horzcat(AeroPanel{:});
                %Determine the AoA distribution along each panel segment
                %
                set(AeroPanel, {'NSPAN'} , num2cell(nSpan));
                set(AeroPanel, {'NCHORD'}, num2cell(nChord));
                set(AeroPanel, {'X1'}    , X1');
                set(AeroPanel, {'X4'}    , X4');
                set(AeroPanel, {'CHORD'} , num2cell(chrd', 1)');
                
            end
            
        end                
        function Dev = checkControlSurfaceOverlap(obj, ControlSurf)
            %checkControlSurfaceOverlap Ensures that there are no control
            %surfaces whose geometry interesects.
            %
            % TODO - This could be done better. Should just get the coords
            % of each coordinate system vertex and then check if one set of
            % coordinates is inside one of the control surfaces.
            
            Dev = [];
            
            if nargin < 2
                ControlSurf = grabControlSurface(obj, 'awi.model.ControlSurface');
            end
            
            if isempty(ControlSurf) %Escape route                
                return
            end
            
            %Identify LE & TE devices by checking their normalised 'x'
            %position and remove any control surfaces that do not belong to
            %either of these sets
            leDev = ControlSurf(arrayfun(@(cs) all(cs.xLE == 0), ControlSurf));
            teDev = ControlSurf(arrayfun(@(cs) all(cs.xTE == 1), ControlSurf));
            Dev   = [leDev ; teDev];
            
            %Check if any of the control surfaces overlap
            %   - There will be an overlap if the start eta position is
            %     less than the end eta position of the previous control
            %     surface.
            %   - Must first ensure they are arranged in order
            leEta = vertcat(leDev.Eta);
            teEta = vertcat(teDev.Eta);
            %
            if ~isempty(leEta)
                [~, ind] = sort(leEta(:, 1));
                leEta = leEta(ind, :);
                leOvrlp = leEta(2 : end, 1) < leEta(1 : end - 1, 2);
            else
                leOvrlp = false;
            end
            if ~isempty(teEta)
                [~, ind] = sort(teEta(:, 1));
                teEta = teEta(ind, :);
                teOvrlp = teEta(2 : end, 1) < teEta(1 : end - 1, 2);
            else
                teOvrlp = false;
            end
            clear ind
            %
            %TODO - Update the error message to tell the user which two
            %control surfaces overlap.
            if or(any(leOvrlp), any(teOvrlp))
                error(['Overlapping coordinate systems found in the ' , ...
                    'LiftingSurface object ''%s''. Unable to convert ', ...
                    'geometry model to analysis model.'], obj.Name);
            end
            
        end          
    end
    
    methods % defining structural layout
        function [Ribs, Stringers] = defineStructuralLayout(obj, Spars)
            %defineStructuralLayout Generates the rib & stringer layout of
            %a generic LiftingSurface torque box based on a cross-section
            %distribution and specified values of rib and stringer pitch
        
            Ribs      = [];
            Stringers = [];
            
            %Parse
            span = obj.Span;
            if isempty(span)
                return
            end
            [rp, sp, Profiles, Spars] = i_parse(obj, Spars);
            
            function [rp, sp, Profiles, Spars] = i_parse(obj, Spars)
                
                rp       = obj.RibPitch;
                sp       = obj.StringerPitch;
                Profiles = obj.CrossSection;
                assert(~isempty(rp), ['The ''RibPitch'' must be specified ', ...
                    'in order to define the structural layout.']);
                assert(~isempty(sp), ['The ''StringerPitch'' must be specified ', ...
                    'in order to define the structural layout.']);
                assert(~isempty(Profiles), ['The LiftingSurface cross-sections ', ...
                    'must be specified in order to define the structural layout.']);
                validateattributes(Spars, {'awi.model.Spar'}, {'nonempty', 'vector'}, ...
                    class(obj), 'Spars');
                %   - Check spars are populated
                sparEta  = get(Spars, {'Eta'});
                sparXLoc = get(Spars, {'XLoc'});
                assert(all(~cellfun(@isempty, sparEta)), ['The ''Eta'' ', ...
                    'property of the Spar objects must be populated in order ', ...
                    'to define the LiftingSurface structural layout.']);
                assert(all(~cellfun(@isempty, sparXLoc)), ['The ''XLoc'' ', ...
                    'property of the Spar objects must be populated in order ', ...
                    'to define the LiftingSurface structural layout.']);
                
                %Order spars from leading edge to trailing edge
                xLoc0 = cellfun(@(x) x(1), sparXLoc);
                [~, ind] = sort(xLoc0, 'ascend');
                Spars = Spars(ind);
                
                %Check that the spars do not cross over at all
                [xLoc_, ~, ~] = interpolateSparXLoc(Spars);
                assert(all(all(diff(xLoc_) > 0)), ['This method assumes ', ...
                    'that none of the spar object cross-over. Update '  , ...
                    'code or implement a bespoke distribution of Rib & ', ...
                    'Spar objects in order to define the structural layout.']);
                
            end
            
            %Define rib locations
            %   - Include control points of aero profile (chord, twist,
            %   etc.) to make sure we capture the planform shape
            etaRib = [0 : rp/span : 1, 1, [Profiles.BeamEta]];
            
            %Get the LE & TE xLoc
            [xLoc, ~, ~, uEta] = interpolateSparXLoc(Spars, obj.AeroEta); %, etaRib);            
            
            %TODO - Extend the use of the Compartment object to allow
            %definition of stringers in the same way as Spars
            
            %Use temp 'Compartment' object to calculate rib cross-sections
            nComp = size(xLoc, 1) - 1;
            Comp  = arrayfun(@(~) awi.model.Compartment, 1 : nComp);
            Crs   = cell(1, nComp);
            for iC = 1 : nComp 
                addCompartment(obj, Comp(iC), ...
                    'EtaLocations', uEta', ...
                    'LocalXStart' , xLoc(iC, :)', ...
                    'LocalXEnd'   , xLoc(iC + 1, :)');
                generateCrossSections(Comp(iC));
                Crs{iC} = Comp(iC).CompartmentCrossSections;
                Crs{iC} = Crs{iC}(ismember([Crs{iC}.BeamEta], uEta));
                removeCompartment(obj, Comp(iC));
            end
            
            %Return the cross-section data
            Ribs = vertcat(Crs{:});            
            
        end
    end
    
    methods % analytic methods  
        function BeamLoc = calculateBeamLoc(obj, eta, chord, xLE)
            %Calculates the normalised position of the beam in the chord
            %direction at the normalised spanwise positions 'eta' using the
            %chord values 'chord'.
            
            %If obj.XLE & obj.XTE have not been defined then default to the
            %interpolation method
            if isempty(xLE)
                BeamLoc = getBPV(obj, 'BeamLoc');
                return
            end
            
            %Logical indexing in case obj.BeamLoc_eta ~= obj.Eta_ (charles added round)
            idx = ismember(round(eta,4), round(obj.BeamLoc_eta,4));
            
            %Shift in global x-direction
            dx  = xLE(idx) + (obj.BeamLoc .* chord(idx));
            
            %Shift in global x-direction at 'Eta_'
            dx_ = awi.model.Component.interp1(obj.BeamLoc_eta, dx, eta);
            BeamLoc = (dx_ - xLE) ./ chord;
            
        end                
        function [MAC, MACxyzLE, MACxyzTE, MACxyzAC] = calculateMAC(obj)
            %calculateMAC Calculate the value and position of the Mean 
            %Aerodynamic Chord (MAC) using equation A-15 (pg 440) from
            %Torenbreek. 
            %
            % The (x,y,z) positions of the MAC are calculated by
            % interpolating the LiftingSurface.
            %
            % The aerodynamic centre of the MAC is assumed to be the
            % quarter-chord.
            
            %Sensible defaults
            MAC      = [];
            MACxyzLE = nan(1, 3);
            MACxyzTE = nan(1, 3);
            MACxyzAC = nan(1, 3);
            
            return
            
            %Grab 'y' positions - Note, not necessarily the coordinates in
            %the y-axes, it depends on the span vector of the lifting
            %surface but the formula should still work.
            y  = obj.SpanVec_;
            c  = obj.Chord_;
            S  = obj.SurfaceArea;
            Si = obj.SegArea;
            tr = obj.TaperRatio;
            
            if isempty(y) || isempty(c) || isempty(S) || isempty(Si) || isempty(tr) %Escape route
                return
            end
            
            %Segment-wise MAC & y-bar
            mac_i  = (2/3) .* ((1 + tr + tr.^2) .* c(1 : end - 1)) ./ (1 + tr);
            ybar_i = obj.Span .* ((1 + 2 .* tr) ./ (3 .* (1 + tr)));
            
            %MAC for the wing
            MAC  = dot(mac_i , abs(Si)) / S;
            ybar = dot(ybar_i, abs(Si)) / S; %Dario's code...
            
            %Calculate y-bar - TODO: Replace this with Dario's code
            %   - Use linear interpolation
            %   - Which segment is the MAC in?
            if numel(c) == 2
                if range(c) == 0    %Special case of untapered wing
                    y_bar = y(1) + 0.5 * (y(end) - y(1));
                else
                    y_bar = interp1(c, y, obj.MAC);
                end
            else
                index = find(MAC < c, 1, 'last');
                c_ = c([index, index + 1]);
                y_ = y([index, index + 1]);
                if range(c) == 0    %Special case of untapered wing
                    y_bar = y(1) + 0.5 * (y(end) - y(1));
                else
                    y_bar = interp1(c_, y_, MAC);
                end
            end
            
            %Find the (x,y,z) position of the beam at the MAC 
            switch obj.SpanVector
                case 'X'
                    error('Fix the implementation of the LiftingSurface for ''SpanVector'' = ''X''.');
                case 'Y'
                    x  = interp1(obj.YData, obj.XData, y_bar);
                    y  = y_bar;
                    z  = interp1(obj.YData, obj.ZData, y_bar);
                    xB = interp1(obj.YData, obj.BeamLoc_i, y);
                case 'Z'
                    x = interp1(obj.ZData, obj.XData, y_bar);
                    y = interp1(obj.ZData, obj.YData, y_bar);
                    z = y_bar;
                    xB = interp1(obj.ZData, obj.BeamLoc_i, y);
            end
            
            %Check for nan
            y(isnan(y)) = 0;            
            
            %Calculate the LE and AC coordinates in the local frame            
            MACxyzLE = [x - (xB .* MAC), y, z];
            MACxyzTE = [x + ((1 - xB) .* MAC), y, z];
            MACxyzAC = [MACxyzLE(1) + (0.25 * MAC), y, z];
            
        end        
        function [cs, eta] = createCoordSysObjects(obj, varargin)
            %createCoordSysObjects Creates the 'awi.model.CoordSys' objects
            %along the span of the LiftingSurface.
            %
            % Ideally we would use the sweep, dihedral & AoA angles to
            % define the orientation of the beam, however, correctly
            % recovering these angles is proving difficult. 
            %
            % For now, just use the coordinate of the trailing edge to
            % define a vector in the XY plane of the beam.
            
            %Allow specification of 'GridVector' & 'GridPlane' values
            p = inputParser;
            addParameter(p, 'GridVector', 'x' , @(x)validatestring(x, {'x', 'y', 'z'}));
            addParameter(p, 'GridPlane' , 'xy', @(x)validatestring(x, {'xy', 'xz', 'yz'}));
            parse(p, varargin{:});
            
            %Create param-value arguments
            args = {'GridVector', p.Results.GridVector, 'GridPlane', p.Results.GridPlane};
            
            %What are the eta positions of the kinks & the chord value at
            %these points?
            chrd = obj.Chord_;
            bLoc = obj.BeamLoc_;
            
            %Offset in global x-direction to the trailing edge
            teOffset = (1 - bLoc) .* chrd;
             
            %Construct the 'G2' offset matrix
            nDat = numel(obj.XData);
            dG2  = [teOffset ; zeros(2, nDat)];
            
            %Use the superclass method
            [cs, eta] = createCoordSysObjects@awi.model.Beam(obj, dG2, args{:});   
            
%             %Orientation angles at the kink locations?       
%             angleP = [obj.Dihedral_ ; obj.AoA_ ; obj.Sweep_]';
%             
%             %Fudge for negative span elements
%             if obj.Span < 0
%                 angleP(:, 1) = angleP(:, 1);
%                 angleP(:, 3) = 180 - angleP(:, 3);
%             end
%             
%             %Fudge for vertically orientated elements
%             if strcmp(obj.SpanVector, 'Z')
%                 aoa = angleP(:, 2);
%                 angleP(:, 2) = angleP(:, 3);
%                 angleP(:, 3) = aoa;
%                 eX = [1 ; 0 ; 0];
%                 eY = [0 ; 1 ; 0];
%                 eZ = [0 ; 0 ; 1];
%                 RI = [eX, eZ, eY];
%             else
%                 RI = eye(3);
%             end
%             
%             %Define 3x3 rotation matrices
%             cs = arrayfun(@(i) awi.model.CoordSys, 1 : nKink, 'Unif', false); 
%             cs = horzcat(cs{:});
%             set(cs, 'ActiveSet', 'pSet');
%             set(cs, 'InitialRMatrix', RI);
%             set(cs, {'AngleP'} , num2cell(angleP, 2));
%             set(cs, {'XOffset'}, num2cell(kinkCoords(:, 1)));
%             set(cs, {'YOffset'}, num2cell(kinkCoords(:, 2)));
%             set(cs, {'ZOffset'}, num2cell(kinkCoords(:, 3)));
%             arrayfun(@(i) build(cs(i)), 1 : nKink);
%             
%             %Store the coordinate systems in the 'Orientation' property of
%             %the underlying 'awi.model.Beam' class 
%             obj.replaceBeamObject(cs, eta);
    
        end            
        function [Crs, eta] = createCrossSectionObjects(obj, varargin)
            %createCrossSectionObjects Generates the cross-section objects
            %that describe the aerofoil profile along the beam and assigns
            %them to the beam.
            %
            % The profiles can be assigned without the model being 'built'
            % but the orientation cannot be defined as it depends on the
            % beam object being fully populated.            
            
            [Crs, ~] = createCrossSectionObjects@awi.model.Beam(obj);
            eta      = obj.AeroEta;
            if isempty(obj.Aerofoil)
                return
            end
                        
            if numel(eta) ~= numel(Crs)
                %Make the 'awi.model.CrossSection' objects and assign the
                %'awi.model.CrossSection' objects to the object
                Crs = arrayfun(@(~) awi.model.CrossSection, 1 : numel(obj.AeroEta));                
                assignBeamObject(obj, Crs, eta, 'replace');                
            end
            CSys = [Crs.Orientation];
            
            %Set up cross-section data
            i_generateAerofoilCoords(obj, eta, Crs);
            
            function i_generateAerofoilCoords(obj, etaQ, CrsSection)
                %i_generateAerofoilCoords Defines the FULL-SCALE
                %coordinates for the aerofoil cross-sections based on the
                %chord distribution and chosen aerofoil profiles.
                
                %Define the normalised coordinates using desired profiles.
                profileNames = obj.Aerofoil;
                Profiles     = arrayfun(@(~) awi.model.CrossSection, 1 : numel(profileNames));
                generateCoordsFromNACA(Profiles, profileNames);
                
                %Interpolate the data to the desired 'eta' positions
                x = vertcat(Profiles.X);
                z = vertcat(Profiles.Z);
                [xU, xL, zU, zL, ~] = interpolateProfileData( ...
                    x, z, obj.Aerofoil_eta, etaQ, Profiles(1).NumPointsX);
                
                %Stash a copy of the coordinate in the 'local' frame for the
                %'CrossSection' data
                xCrs = [xU, fliplr(xL)];
                zCrs = [zU, fliplr(zL)];
                
                %Grab planform parameters
                planformEta = getEta(obj , 'Planform Property', 'AxisFlag', obj.SpanVector);
                chrd        = getBPV(obj, 'Chord');
                
                %Interpolate planform parameters to eta positions of
                %profiles
                chord   = interp1(planformEta, chrd, etaQ);
                xLE     = awi.model.LiftingSurface.interp1(planformEta, obj.XLE, etaQ);
                beamLoc = calculateBeamLoc(obj, etaQ, chord, xLE);
                
                %Set up for implicit array expansion
                chord   = chord';
                beamLoc = beamLoc';
                
                %Scale the (x, z) coordinates and shift by the beam offset
                xCrs = (xCrs .* chord) - (beamLoc .* chord);
                zCrs = zCrs .* chord;
                
                set(CrsSection, {'X'}, num2cell(xCrs, 2), {'Z'}, num2cell(zCrs, 2));
                
                %Set profile names where possible
                idx_ = ismember(etaQ, obj.Aerofoil_eta);
                set(CrsSection(idx_), {'Name'}, {Profiles.Name}');
                
            end
            
            %% Assign the coordinate system transforms as well!
            %   - Account for AoA and dihedral but no sweep
            
            %Define the 'Origin' of the coordinate system as the beam
            %locations
            [xd, yd, zd] = xyzdata(obj);
            if isempty(xd)
                return
            end
            beam = obj.AbsPosition' + [xd ; yd ; zd];
            set(CSys, {'Origin'}, num2cell(beam', 2));
                        
            %Vector pointing along the beam axis
            v     = diff(beam, [], 2);
            v_mod = sqrt(sum(v.^2, 1));
            v     = v ./ repmat(v_mod, [3, 1]);
            
            %Assume end coordinate system has same orientation
            v = [v, v(:, end)];
            
            %Remove effects of sweep - based on the 'SpanVector' value
            switch obj.SpanVector
                case 'X'
                    error('Update code for ''SpanVector'' = ''X''.');
                case 'Y'
                    v(1, :) = 0; 
                    ind = 2;
                case 'Z'
                    v(1, :) = 0; 
                    ind = 3;
            end
            
            %Account for flipped rotation matrix for portside objects
            idx = v(ind, :) < 0;
            v(ind, idx) = -v(ind, idx);
            
            %Define 'gSet' quantities
            G1 = (beam + v)';
            G2 = obj.TE';
            gridVec   = 'x';
            gridPlane = 'xy';
            
            %Assign 'gSet' quantities and build the coordinate systems
            set(CSys, 'ActiveSet', 'gSet');
            set(CSys, 'GridVector', gridVec);
            set(CSys, 'GridPlane' , gridPlane);
            set(CSys, {'G1'}      , num2cell(G1, 2));
            set(CSys, {'G2'}      , num2cell(G2, 2));            
            arrayfun(@build, CSys);
            
            %% Account for folding wing tips
            
            %Grab the coordinate system of the at the hinge point
            [HingeCoordSys, flareRotMatrix] = grabHingeCoordSys(obj);
            
            %Update orientation of the cross-section at this point
            %   - TODO: Account for 'shortening' of the chord
            HingeCoordSys.ActiveSet = 'rSet';
            HingeCoordSys.RMatrix   = flareRotMatrix;
            build(HingeCoordSys);
            
        end        
    end
    
    methods % shortcut for grabbing control surfaces 
        function ControlSurfaceObject = grabControlSurface(obj, type_or_class)
            %grabControlSurface Extracts any child control surface objects
            %corresponding to the type or class specified in
            %'type_or_class'.
            
            %Start with an empty
            ControlSurfaceObject = [];
            
            %Escape route
            if isempty(obj.Children)
                return
            end
            
            %Token for generic control surface
            genTok = {'generic' ; 'awi.model.ControlSurface'};
            
            %Construct valid names
            cls = cellfun(@(fn) func2str(fn), obj.ControlSurfaceSubTypes(:, 1), 'Unif', false);
            validStrCell = [cls, obj.ControlSurfaceSubTypes(:, [2, 3])];
            validString  = [genTok ; validStrCell(:)];
            
            %Validate input
            validatestring(type_or_class, validString, class(obj), 'the requested control surface type');
            
            %Search for the token
            switch type_or_class
                case genTok
                    tok = 'awi.model.ControlSurface';
                otherwise
                    idx = arrayfun(@(i) any(strcmpi(validStrCell(i, :), type_or_class)), 1 : size(validStrCell, 1))';                    
                    tok = cls{idx};                    
            end
            
            %Children of this liftng surface
            ch = obj.Children;
            
            %Only retain collector objects otherwise we will get control
            %surfaces that belong to child lifting surfaces
            idx = arrayfun(@(x) isa(x, 'mvc.model.Collector'), ch);
            
            %Grab the control surfaces corresponding to the specified type           
            allCh = flatlist(obj.Children(idx));            
            idx   = arrayfun(@(a) isa(a, tok), allCh);
            ControlSurfaceObject = allCh(idx);
            
            %Force empty matrix if we cannot find the specified object
            if isempty(ControlSurfaceObject)
                ControlSurfaceObject = [];
            end
            
        end        
    end
    
    methods % generating camber data
        function AerofoilData = getAerofoilData(obj, AerofoilCoords)
            %getAerofoilData Structure containing a breakdown of each
            %aerofoil section. Details include:
            %   - Upper &  lower surface coordinates
            %   - Mean camber line
            %   - Camber angle
            
            %Grab data
            x = AerofoilCoords.X;
            y = AerofoilCoords.Y;
            z = AerofoilCoords.Z;
            
            %Check that the aerofoil coordinates are situated in a uniform,
            %normalised grid
            if any(any(x > 1)) || any(any(x < 0)) || any(any(y > 1)) || any(any(y < 0))
                ME = MException('matlab:awi:badAerofoilCoords'        , ...
                    ['The aerofoil coordinates in the structure '     , ...
                    '''AerofoilCoords'' for the object ''%s'' must '  , ...
                    'be normalised in the chorwise and spanwise '     , ...
                    'direction.'], obj.NameAndType);
                throw(ME);
            end
            
            %Check that our assumption is valid.
            if size(x, 1) ~= 2*obj.NumAerofoilPointsX
                error(['Assumption that we can index the aerofoil ', ...
                    'coordinates to extract the upper and lower '  , ...
                    'surfaces is invalid. Add new code']);
            end
            
            %Seperate into upper and lower coordinates
            %   - The assumtpion here is that the coordinates have already
            %     been preprocessed and therefore they have the same number
            %     of coordinates on the upper and lower surfaces.
            %   - A more generic check may be required in future. An
            %     improvement would be to find the leading edge and
            %     trailing edge of the aerofoil and then split at that
            %     point.
            AerofoilData.XUpper = x(1 : obj.NumAerofoilPointsX, :);
            AerofoilData.XLower = x(obj.NumAerofoilPointsX + 1 : end, :);
            AerofoilData.ZUpper = z(1 : obj.NumAerofoilPointsX, :);
            AerofoilData.ZLower = z(obj.NumAerofoilPointsX + 1 : end, :);
            
            %Find the mean camber line - TODO This is not corrrect, the
            %camber line is the line between the locus of points that are
            %equidistant from the upper and lower surfaces.
%             AerofoilData.xCamber = AerofoilData.XUpper;
%             AerofoilData.zCamber = mean(cat(3, AerofoilData.ZUpper, AerofoilData.ZLower), 3);
            AerofoilData.xCamber = [];
            AerofoilData.zCamber = [];
            
            %Calculate the camber angle
            AerofoilData.xCamberAngle = [];
            AerofoilData.CamberAngle  = [];
            
        end    
    end
    
    methods % folding wing tips
        function FWT = insertWingFold(obj, varargin)
            %insertWingFold Inserts a wing-fold at a specified position
            %along the span with a specified flare angle.
            %
            % Parameters:
            %   - 'EtaFold'   : Spanwise position of the fold (normalised)
            %   - 'FlareAngle': Flare angle of the folding wing tip
            %
            % Detailed Description:
            %   - Splits the current LiftingSurface at the desired fold
            %     positions and adds a new LiftingSurface object to the tip
            %     of the wing. All beam/planform properties are split at
            %     the wing fold position and reassigned to the inboard/FWT
            %     objects to account for the new 'Span' of the wing.
            %
            % Assumptions:
            %   - Assume the same material is used for the inboard and
            %     outboard wing. This will not be valid if the
            %     LiftingSurface has multiple Material objects along the
            %     span.
            %
            % TODO - Resolve issues with interpolating the cross-section
            % data using 'getBPV'
            
            FWT = []; % The Folding-Wing-Tip (FWT) as a new LiftingSurface
            
            if (numel(obj) > 1)
                warning(['Method ''insertWingFold'' is not valid for ', ...
                    'object arrays. Returning to invoking function.']);
                return
            end
            
            p = inputParser;
            addParameter(p, 'EtaFold', 0.75, @(x)validateattributes(x, ...
                {'numeric'}, {'scalar', 'positive', '<=', 1, 'nonnan', ...
                'finite', 'real'}, class(obj), 'EtaFold'));
            addParameter(p, 'FlareAngle', 0, @(x)validateattributes(x, ...
                {'numeric'}, {'scalar', 'nonnan', 'finite', 'real'}, ...
                class(obj), 'FlareAngle'));
            addParameter(p, 'FoldAngle', 0, @(x)validateattributes(x, ...
                {'numeric'}, {'scalar', 'nonnan', 'finite', 'real'}, ...
                class(obj), 'FoldAngle'));
            parse(p, varargin{:});
            
            FWT = awi.model.LiftingSurface('Name', [obj.Name, '_FWT']);
            FWT.SpanVector       = obj.SpanVector;
            FWT.IsFoldingWingTip = true;
            
            %Define eta positions in 'Span' and 'R' axes for splitting data
            % - We have to make a distinction between beam properties which
            % are defined along the span direction and those which are
            % defined along the piecwise straight sections of the beam
            etaFold = p.Results.EtaFold;
            switch obj.SpanVector
                case 'Y'
                    s_data = obj.YData;
                case 'Z'
                    s_data = obj.ZData;
            end
            r_data = obj.RData;            
            span   = obj.Span;
            length = r_data(end);
            etaS   = s_data/span;
            etaR   = r_data / length;
            etaR   = interp1(etaS, etaR, etaFold);
            
            %Split all beam properties at the desired fold position
            bpAxis = {obj.BeamProperties.AxisFlag};
            indexR = find(and(ismember(bpAxis, 'R'), ismember({obj.BeamProperties.Type}, 'FE BeamProp')));
            indexS = find(ismember(bpAxis, obj.SpanVector));
            [obj, FWT] = i_assignBeamPropData(obj, FWT, etaFold, [indexR,indexS]);
            
            function [obj, FWT] = i_assignBeamPropData(obj, FWT, etaF, index)
                %i_assignBeamPropData Reassigns the beam property data to
                %account for the new LS length that results from adding a
                %FWT.
                %
                % TODO - We should be able to link directly to the
                % BeamProperties and the pre-post listeners should
                % correctly assign the data but that doesn't seem to be
                % working at the moment. Instead we use dynamic property
                % referencing...
                
                for ii = index
                    
                    nam_str = obj.BeamProperties(ii).Quantity;
                    eta_str = [nam_str, '_eta'];
                     
                    %Add a new data point at the desired fold position
%                     eta   = obj.BeamProperties(ii).Distribution;
%                     val   = obj.BeamProperties(ii).Value;
                    eta = obj.(eta_str);
                    val = obj.(nam_str);
                    bpVal = getBPV(obj.BeamProperties(ii), nam_str, etaF);
                    [eta_, ia] = unique([eta, etaF]);
                    val_ = [val, bpVal];
                    val_ = val_(ia);
                    
                    %Split into inboard and outboard sections
                    idx_i = eta_ <= etaF;
                    idx_o = eta_ >= etaF;
                    
                    eta_i = eta_(idx_i);
                    obj.(nam_str) = val_(idx_i);
                    obj.(eta_str) = eta_i / eta_i(end);
%                     obj.BeamProperties(ii).Distribution = eta_i / eta_i(end);
%                     obj.BeamProperties(ii).Value        = val_(idx_i);
                    
                    eta_o = eta_(idx_o);
                    FWT.(nam_str) = val_(idx_o);
                    FWT.(eta_str) = (eta_o - eta_o(1)) / (eta_o(end)- eta_o(1));
%                     FWT.BeamProperties(ii).Distribution = (eta_o - eta_o(1)) / (eta_o(end)- eta_o(1));
%                     FWT.BeamProperties(ii).Value        = val_(idx_i);
                    
                end
                
            end
            
            %TODO - Add call to 'interpolateCrossSection' here?
            obj.FlareAngle = p.Results.FlareAngle;
            FWT.FlareAngle = p.Results.FlareAngle;
            obj.FoldAngle  = p.Results.FoldAngle; 

            %Assume the same material
            FWT.Material     = obj.Material;
            FWT.Material_eta = [0, 1];
            
            %Update the geometry
            FWT.Span    = (1 - etaFold) * span;
            obj.Span    = etaFold * span;
            
            chord=obj.Chord;
            FWT.RootChord=chord(end);

            %Add the FWT to the end of the inboard wing
            FWT.SOffset = 1;
            obj.add(FWT);            
            %Propogate changes
            FWT.ActiveSet = obj.ActiveSet;
            build(obj);

        end
        function [HingeCoordSys, flareRotMatrix] = grabHingeCoordSys(obj)
            %grabHingeCoordSys Selects the hinge CoordSystem object and
            %rotates it about a vector which is normal to the sweep plane
            %by the desired flare angle.
            
            Crs = obj.CrossSection;
            
            HingeCoordSys  = [];
            flareRotMatrix = [];
            
            if obj.HasFoldingWingTip
                HingeCoordSys = Crs(ismember([Crs.BeamEta], 1)).Orientation;
            elseif obj.IsFoldingWingTip
                HingeCoordSys = Crs(ismember([Crs.BeamEta], 0)).Orientation;
            else
                return
            end
            
            %Rotate the coordinate system (about the axis which is normal 
            %to the sweep plane) by the flare angle
            switch obj.SpanVector 
                case 'Y'
                    rot_vec = [0, 0, 1];
                case 'Z'
                    rot_vec = [0, 1, 0];
            end
            rot_angle = -obj.FlareAngle; %right-hand corkscrew rule
            rotMatrix = awi.model.CoordSys.rodriguezRotation(rot_vec, rot_angle);

%             if obj.FlareAngle > 0
%                 rotation_dihedral=awi.model.CoordSys.rodriguezRotation([0,1,0], obj.Dihedral(1));
%             else
%                 rotation_dihedral=awi.model.CoordSys.rodriguezRotation([0,1,0], -obj.Dihedral(1));
%             end


            rotMatrix_inv = awi.model.CoordSys.rodriguezRotation(rot_vec, -rot_angle);
            foldMatrix = awi.model.CoordSys.rodriguezRotation([1,0,0], obj.FoldAngle);
%             flareRotMatrix = HingeCoordSys.RMatrix' * rotMatrix * foldMatrix * rotMatrix_inv;
            flareRotMatrix = rotMatrix*HingeCoordSys.RMatrix';
        end
    end
    
    methods (Static) %angle2Beam, shoelaceFormula, getAeroPanelSize
        
        function gridVector = angle2Beam(spanVec, angle)
            %angle2Beam Calculates a vector of grid points based on the
            %position of the points along a given span vector and the angle
            %of the vector to that span vector.
            %
            % The angle must be in units of degrees [deg].
            %
            % It is assumed that the coords start at [0,0,0].
            
            gridVector = cumsum([0, diff(spanVec) .* tand(angle(1:end-1))]);
            
        end        
        
        function A = shoelaceFormula(coordinates, varargin)
            %shoelaceFormula Uses the shoelace formula to calculate the
            %area of a 2D polgyon described by coordinates. More
            %information about the shoelace formula can be found <a
            %href="https://en.wikipedia.org/wiki/Shoelace_formula">here</a>.
            %
            % A = shoelaceFormula(coords, 'XY') calculates the
            % area of the polygon described by the coordinates 'coords' in
            % the plane 'XY'.
            %
            % The number of columns in 'coords' cannot exceed 3. It is
            % assumed that the order of the columns in 'coords' is [x,y,z].
            %
            
            %parse inputs
            p = inputParser;
            addRequired(p, 'coordinates', @(x)validateattributes(x, ...
                {'numeric'}, {'nonnan', 'nonempty', 'finite', '2d'}));
            addOptional(p, 'Plane', 'XY', @(x)any(validatestring(x, ...
                {'XY', 'XZ', 'YZ'})));
            parse(p, coordinates, varargin{:});
            
            % check size of 'coordinates'
            nCol = size(coordinates, 2);
            if nCol > 3 || nCol < 2
                ME = MException( ...
                    'MATLAB:uob:LiftingSurface:BadPlanformCoordSize' , ...
                    ['The number of columns in ''%s'' cannot exceed ', ...
                    '3 and must be greater than 1. The current '     , ...
                    'number of columns is %i.'], 'coordinates', nCol);
                throwAsCaller(ME);
            end
            
            % grab sets of coordinates
            if nCol == 3
                % if (x,y,z) coordinates have been provided decide which
                % set of coordinates to use ...
                switch p.Results.Plane
                    case 'XY'
                        r = coordinates(:, 1);
                        s = coordinates(:, 2);
                    case 'XZ'
                        r = coordinates(:, 1);
                        s = coordinates(:, 3);
                    case 'YZ'
                        r = coordinates(:, 2);
                        s = coordinates(:, 3);
                end
            else
                % if only 2 columns have been provided then grab data
                r = coordinates(:, 1);
                s = coordinates(:, 2);
            end
            
            %TODO - Check that the coordinates are in
            %clockwise/anti-clockwise order!
            
            % calculate the area using the shoelace formula [1].
            A = 0.5 .*abs( ...
                sum(r(1:end-1, :) .* s(2 : end, :), 1) - ...
                sum(r(2:end, :)   .* s(1:end-1, :), 1) + ...
                r(end, :) .* s(1, :)                   - ...
                r(1, :)   .* s(end, :));
        end
        
        function [nSpan, nChord] = getAeroPanelSize(obj, span, chord)
            %getAeroPanelSize Returns the number of spanwise and chordwise
            %panels in each aerodynamic panel defined by 'span' & 'chord'.
            
            %Calculate aero panel size
            %   - 'panelLength' is dimension in the global X direction
            %   - 'panelWidth' is dimension in the global spanwise axis
            if isempty(obj.AeroPanelLength)
                panelLength = chord ./ obj.NumAeroPanel;
            else
                panelLength = repmat(obj.AeroPanelLength, size(span));
            end
            panelWidth = panelLength .* obj.AeroPanelAR;
            
            %Determine number of panels in the chordwise and spanwise direction
            nSpan   = ceil(span ./ panelWidth);
            nChord  = ceil(chord ./ panelLength);
            
        end
        
    end

    methods (Static) %Default builds
        function LSObj = makeHodgesWing 
            %makeHodgesWing Makes a lifting surface object representing the
            %starboard-wing of the Hodges-Patil HALE UAV.
            % 
            % Properties:
            %   - Semi-Span : 16     m
            %   - Chord     : 1      m
            %   - NSM       : 0.75   kg/m
            %   - NSI       : 0.1    kgm
            %   - EI_xx     : 2x10^6 Nm^2
            %   - EI_zz     : 4x10^6 Nm^2
            %   - GJ        : 1x10^4 Nm^2
            %   - Beam Loc  : 50% chord 
            %
            % Notes:
            %   - Neutral & shear axis are coincident, CoM lies on beam.
            %   - Specify a large cross-sectional area so the beam if
            %     effectively stiff in the axial DOF.
            %
            % For more information see: https://doi.org/10.2514/2.2738
          
            LSObj = awi.model.LiftingSurface;
            LSObj.Name = 'HodgesWing';
            
            %As awi.model.Beam uses I11/I22/J not EI/GJ we must specify a
            %notional material which yields I/J from EI/GJ.
            E  = 76e9; %[N/m^2], typical YM of aluminium 
            nu = 0.333;
            
            %Make the material
            Mat = awi.model.Material;
            Mat.E  = E;
            Mat.Nu = nu;
            Mat.G  = E / (2 * (1 + nu));
            
            %Use the 'Parametric' set for building the wing
            LSObj.ActiveSet = 'pSet';
            
            %Geometry
            LSObj.SpanVector   = 'Y'; %standard Nastran coordinate system
            LSObj.Span         = 16;
            LSObj.Chord        = [1, 1];
            LSObj.Chord_eta    =  [0, 1];
            LSObj.Aerofoil     = {'NACA0012', 'NACA0012'};
            LSObj.Aerofoil_eta = [0, 1];
            LSObj.BeamLoc      = 0.5;
            
            %Beam Properties
            mybox=awi.model.BoxBeam;
            mybox.BoxType='SymmetricBox';
            mybox.Height=0.2;
            mybox.Width=1;
            mybox.CoverThickness=0.001;
            mybox.SparThickness=0.002;
            getGeometricProps(mybox)
            LSObj.BoxBeam = mybox;
            LSObj.A   = mybox.Abb;    %dummy value, required by Nastran to run solution
            LSObj.I11 = mybox.Ixx;
            LSObj.I22 = mybox.Izz;
            LSObj.J   = mybox.Jbb;
            LSObj.NSM = mybox.NSM;
            LSObj.NSI = mybox.NSI;
            
            %Material 
            LSObj.Material_eta = [0, 1];
            LSObj.Material     = [Mat, Mat];
            
            %Build the object in order to populate remaining properties
            build(LSObj);
        end
        function LSObj = makeA320Wing
            %makeA320Wing Makes a lifting surface object which is
            %representative (but not an exact match) of the starboard-wing
            %of the Airbus A320.
            %
            % Properties:
            %
            % Notes:
            %   - 
            %
            % For more information see: https://en.wikipedia.org/wiki/Airbus_A320_family#Specifications
            
            LSObj = awi.model.LiftingSurface;
            LSObj.Name = 'A320Wing';
            
            %Use the Leading/Trailing edge sweep to define the planform
            LSObj.ActiveSet = 'sSet';
            
            %Wing dimensions
            LSObj.SpanVector  = 'Y';
            LSObj.Span        = 35.8/2;
            LSObj.LESweep     = [25, 25];
            LSObj.LESweep_eta = [0, 1];
            LSObj.TESweep     = [0, 10, 10];
            LSObj.TESweep_eta = [0, 0.3, 1];
            LSObj.RootChord   = 5;
            
            %Make sure the beam is at the midchord
            all_eta           = LSObj.Eta_;
            LSObj.BeamLoc     = repmat(0.5, size(all_eta));
            LSObj.BeamLoc_eta = all_eta;
            
            %Spars - 2 x spars, 1 @ 15% chord, 1 @ 65% chord
            FrontSpar = awi.model.Spar;
            FrontSpar.XLoc = [0.15, 0.15];
            FrontSpar.Eta  = [0   , 1];
            RearSpar = awi.model.Spar;
            RearSpar.XLoc = [0.65, 0.65];
            RearSpar.Eta  = [0   , 1];
            LSObj.add([FrontSpar, RearSpar]);
            
            %Define internal layout
            LSObj.RibPitch      = 0.65;
            LSObj.StringerPitch = 0.15;

            %Make the material-charles---------
            E  = 76e9; %[N/m^2], typical YM of aluminium 
            nu = 0.333;
            Mat = awi.model.Material;
            Mat.E  = E;
            Mat.Nu = nu;
            Mat.G  = E / (2 * (1 + nu));
            
            LSObj.Material_eta = [0, 1];
            LSObj.Material     = [Mat, Mat];
            
            %Beam property definition - charles
            mybox=awi.model.BoxBeam;
           
            mybox.BoxType='SymmetricBox';
            mybox.Height=0.2;
            mybox.Width=2.5;
            mybox.CoverThickness=0.001;
            mybox.SparThickness=0.002;
            getGeometricProps(mybox)
            LSObj.BoxBeam = mybox;
            
            LSObj.A   = mybox.Abb;    
            LSObj.I11 = mybox.Ixx;
            LSObj.I22 = mybox.Izz;
            LSObj.J   = mybox.Jbb;
            LSObj.NSM = mybox.NSM;
            LSObj.NSI = mybox.NSI;
            
            % Aeropanel - charles
           
            LSObj.NumAeroPanel=30;
            
            %Control surfaces - flaps
            myflap=awi.model.Flap;

            myflap.Eta=[0.1, 0.2];
            myflap.xLE=[0.8,0.8];
            myflap.xTE=[1,1];
            myflap.Max_def=0.1;
            myflap.Max_rate=0.1;
            myflap.HingeLine='LE';
            myflap.Label='flap1';
            build(myflap)
            LSObj.add(myflap);
            
            %Build the object in order to populate remaining properties
            build(LSObj);            
        end
        function LSObj = makeStraightWingWithFold(eta_fold, flare_angle)
            %makeStraightWingWithFold Makes a straight, unswept, untapered
            %wing with a fold.
            %
            % Inputs:
            %   - eta_fold   : Position of the fold along the span
            %   - flare_angle: Angle of the FWT flare
            %
            % Notes:
            %   - Uses the Hodges Wing planform but inserts a wing fold at
            %     desired location with a specified flare angle.
            
            %Fold parameters
            if nargin < 2
                eta_fold    = 0.75; %[-], normalised position of the fold along the semi-span
                flare_angle = 0;    %[-], angle of wing-fold w.r.t global-X
            end
            
            %Make the objects
            LSObj = awi.model.LiftingSurface.makeHodgesWing;
            FWT   = insertWingFold(LS, ...
                'EtaFold'   , eta_fold, ...
                'FlareAngle', flare_angle); %#ok<NASGU>
           
        end
    end
    
end
