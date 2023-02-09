classdef (ConstructOnLoad) ControlSurface < awi.model.Component
    %ControlSurface Defines a generic control surface located on a
    %'LiftingSurface' object.
    
    %Parameters
    properties (AbortSet, SetObservable)
        %Non-dimensional spanwise position of the start and end of the
        %control surface with respect to its parent.
        Eta
        %Normalised chordwise position of the leading edge of the control
        %surface with respect to the leading edge of its parent.
        xLE
        %Normalised chordwise position of the trailing edge of the control
        %surface with respect to the leading edge of its parent.
        xTE
        %Maximum permitted deflection angle of the control surface.
        Max_def
        %Maximum permitted deflection rate of the control surface.
        Max_rate
        %Defines the hinge line of the control surface
        HingeLine = 'LE';
        %Defines the change in lift coefficient per unit deflection of the
        %control surface
        DeltaCL
        %Defines the change in moment coefficient per unit deflection of 
        %the control surface
        DeltaCM
    end
    
    %Appearance of the patch
    properties
        %And face colour and transparency for patch
        FaceColor = 'r';
        FaceAlpha = 1;
    end
    
    %Name of the control surface (short-hand)
    properties
        %Shorthand name of the control surface (For MSC.Nastran representation)
        Label
    end
    
    %Coordinate System - Add this as a child once we have sorted the hidden
    %tree object functionality
    properties (Dependent)
        %3x3 rotation matrix describing the coordinate system of the
        %control system hinge
        HingeCoordSys
    end
    
    %Shadow properties
    properties (Dependent)
        %Corner coordinates of the aerodynamic surface in the global
        %coordinate system
        Coords
        %Absolute position of the parent object
        ParentAbsPos
        %Parent of the ControlSurface
        ParentBeam
    end
    
    methods % set / get
        function val = get.ParentBeam(obj)
            
            val = obj.Parent;
            
            while ~isempty(val) && ~isa(val,'awi.model.Beam')
                
                val = val.Parent;
                
            end
            
        end
        function set.Eta(obj, val)            %set.Eta           
            %set.Eta Set method for the property 'Eta'.
            %
            % Rules :
            %   - 'Eta' must be of type double with the following
            %     attributes: row, 2 elements, nonnan, finite, real,
            %     nonnegative, <= 1.
            
            validateattributes(val, {'double'}, {'row', 'numel', 2, ...
                'finite', 'real', 'nonnan', 'nonnegative', '<=', 1}, ...
                class(obj), 'Eta');
            obj.Eta = val;
        end        
        function set.xLE(obj, val)            %set.xLE           
            %set.xLE Set method for the property 'xLE'.
            %
            % Rules :
            %   - 'xLE' must be of type double with the following
            %     attributes: row, 2 elements, nonnan, finite, real,
            %     nonnegative, <= 1.
            
            validateattributes(val, {'double'}, {'row', 'numel', 2, ...
                'finite', 'real', 'nonnan', 'nonnegative', '<=', 1}, ...
                class(obj), 'xLE');
            obj.xLE = val;
        end        
        function set.xTE(obj, val)            %set.xTE           
            %set.xTE Set method for the property 'xTE'.
            %
            % Rules :
            %   - 'TE' must be of type double with the following
            %     attributes: row, 2 elements, nonnan, finite, real,
            %     nonnegative, <= 1.
            
            validateattributes(val, {'double'}, {'row', 'numel', 2, ...
                'finite', 'real', 'nonnan', 'nonnegative', '<=', 1}, ...
                class(obj), 'xTE');
            obj.xTE = val;
        end        
        function set.Max_def(obj, val)        %set.Max_def       
            %set.Max_def Set method for the property 'Max_def'.
            %
            % Rules :
            %   - 'Max_def' must be of type double with the following
            %     attributes: scalar, nonnan, finite, real.
            
            validateattributes(val, {'double'}, {'scalar', 'finite', ...
                'real', 'nonnan'}, class(obj), 'LE');
            obj.Max_def = val;
        end    
        function set.Max_rate(obj, val)       %set.Max_rate      
            %set.Max_rate Set method for the property 'Max_rate'.
            %
            % Rules :
            %   - 'Max_rate' must be of type double with the following
            %     attributes: scalar, nonnan, finite, real.
            
            validateattributes(val, {'double'}, {'scalar', 'finite', ...
                'real', 'nonnan'}, class(obj), 'LE');
            obj.Max_rate = val;
        end   
        function set.HingeLine(obj, val)      %set.HingeLine     
            %set.HingeLine Set method for the property 'HingeLine'.
            %
            %   - 'obj.HingeLine' must be one of the valid strings 'LE' or
            %   'TE'.
            
            validatestring(val, {'LE', 'TE'}, class(obj), 'HingeLine');
            obj.HingeLine = val;
        end
        function set.Label(obj, val)          %set.Label         
            %set.Label Set method for the property 'Label'.
            %
            % Label must be a character vector of no-more than 8 characters
                        
            if isempty(val) %Force an empty matrix
                obj.Label = [];
                return
            end            
            validateattributes(val, {'char'}, {'row', 'nonempty'}, ...
                class(obj), 'Label');
            assert(numel(val) < 9, ['The control surface label must ', ...
                'be no more than 8 characters in width.']);
            obj.Label = val;
        end
        function val = get.HingeCoordSys(obj) %get.HingeCoordSys 
            %get.HingeCoordSys Get method for the dependent property
            %'HingeCoordSys'.
            %
            %   - 'HingeCoordSys' is an object of type 'awi.model.CoordSys'
            %   - The 3x3 rotation matrix is defined using the 'gSet'
            %     parameter set of the 'awi.model.CoordSys' class.
            %   - The first grid point (G1) defines the local y-axis of the
            %     coordinate system which must point along the hinge line 
            %     of the control surface. This is in keeping with the 
            %     MSC.Nastran convention.
            %   - The second grid point (G2) defines the local XY-plane of
            %     the control surface. The local z-axis is then calculated 
            %     by the 'gSet' build method in the 'awi.model.CoordSys'
            %     class. 
                    
            %Get coordinates of the corners of the control surface in the
            %global coordinate system.
            ControlSurfCoords = obj.Coords;
            
            %Initialise object
            val = awi.model.CoordSys;
            val.ActiveSet  = 'gSet';
            val.GridVector = 'y';
            val.GridPlane  = 'xy';
            
            %Escape route for when 'get' is called on empty object
            if isempty(ControlSurfCoords) || isempty(ControlSurfCoords.(obj.HingeLine).X)
                return
            end
            
            %Set origin of Coordinate System as the inboard corner of the
            %control surface LE
            %   - TO DO : This is not valid for spoilers!!
            O = [ ... 
                ControlSurfCoords.(obj.HingeLine).X(1) , ...
                ControlSurfCoords.(obj.HingeLine).Y(1) , ...
                ControlSurfCoords.(obj.HingeLine).Z(1)];
            val.Origin = O;
            
            %The first grid defines the y-axis of the coordinate system
            val.G1 = [ ... 
                ControlSurfCoords.(obj.HingeLine).X(2) , ...
                ControlSurfCoords.(obj.HingeLine).Y(2) , ...
                ControlSurfCoords.(obj.HingeLine).Z(2)];
            
            %The second grid defines a point in the local XY-plane
            switch obj.HingeLine
                case 'LE'
                    val.G2 = [ ...
                        ControlSurfCoords.TE.X(2) , ...
                        ControlSurfCoords.TE.Y(2) , ...
                        ControlSurfCoords.TE.Z(2)];
                case 'TE'
                    val.G2 = [ ...
                        ControlSurfCoords.LE.X(2) , ...
                        ControlSurfCoords.LE.Y(2) , ...
                        ControlSurfCoords.LE.Z(2)];
            end
            
            %Build the 'awi.model.CoordSys' object using the 'gSet'.
            build(val);

        end            
        function val = get.Coords(obj)        %get.Coords        
            %get.Coords Get method for the property 'Coords'.
            %
            % 'Coords' is a structure containing the x, y, z coordinates of
            % the corners of the control surface. The coordinates are in 
            % global coordinate system and are NOT relative to the origin
            % of the parent LiftingSurface object.
            
            temp = struct('X', [], 'Y', [], 'Z', []);
            val  = struct('LE', temp, 'TE', temp);
            
            %Check if parent is of type 'LiftingSurface'
            ph = obj.Parent;
            if isempty(ph)
                return
            elseif isa(ph, 'mvc.model.Collector')
                ph = ph.Parent;
            end
            
            %Position of origin relative to parent
            O = obj.ParentAbsPos;
            
            %Grab properties of parent LiftingSurface object                   
            planformCoords = ph.PanelCoords;            
            
            %Calculate the chord
            chord = planformCoords.TE.X - planformCoords.LE.X;
            
            %Nothing ?
            if isempty(chord)
                
                %Bail
                val = [];
                return;
                
            end
            
            %Calculate the LE coordinates and chord values at 'obj.Eta' 
            switch ph.SpanVector
                case 'X'
                    sVec = obj.Eta .* ph.Span + O(1);
                    x_LE = obj.xLE;
                    x_TE = obj.xTE;
                    x = sVec;
                    y = interp1(planformCoords.LE.X, planformCoords.LE.Y, sVec);
                    z = interp1(planformCoords.LE.X, planformCoords.LE.Z, sVec);
                case 'Y'
                    %Grab the spanwise coordinates in the global-CS
                    sVec = obj.Eta .* ph.Span + O(2);
                    x_LE = obj.xLE;
                    x_TE = obj.xTE;
                    %Check if the control surface is defined over multiple
                    %planform bays
                    idx = and(abs(planformCoords.LE.Y) > abs(sVec(1)),  abs(planformCoords.LE.Y) < abs(sVec(2)));
                    if any(idx)         
                        %Determine sort method
                        del = sVec(end) - sVec(1);
                        if del < 0
                            order = 'descend';
                        else
                            order = 'ascend';
                        end
                        %Add new points and sort
                        sVec = sort([sVec, planformCoords.LE.Y(idx)], order);
                        %Interpolate the 'xLE' and 'xTE' values assuming a
                        %linear variation
                        x_LE = interp1(obj.Eta, obj.xLE, sVec ./ ph.Span, 'linear', 'extrap');
                        x_TE = interp1(obj.Eta, obj.xTE, sVec ./ ph.Span, 'linear', 'extrap');
                    end
                    %Inteprolate the values at the control surface
                    x = interp1(planformCoords.LE.Y, planformCoords.LE.X, sVec);
                    y = sVec;
                    z = interp1(planformCoords.LE.Y, planformCoords.LE.Z, sVec);
                    c = interp1(planformCoords.LE.Y, chord, y);
                    
                case 'Z'
                    sVec = obj.Eta .* ph.Span + O(3);
                    x_LE = obj.xLE;
                    x_TE = obj.xTE;
                    x = interp1(planformCoords.LE.Z, planformCoords.LE.X, sVec);
                    y = interp1(planformCoords.LE.Z, planformCoords.LE.Y, sVec);
                    z = sVec;
                    c = interp1(planformCoords.LE.Z, chord, z);
            end

            %Define absolute coordinates of the LE & TE corners of the
            %control surface
            val.LE.X = x + (x_LE .* c);
            val.LE.Y = y;
            val.LE.Z = z;
            val.TE.X = x + (x_TE .* c);
            val.TE.Y = y;
            val.TE.Z = z;
            
        end   
        function val = get.ParentAbsPos(obj)  %get.ParentAbsPos  
            %get.ParentAbsPos Get method for the dependent property
            %'ParentAbsPos'
            %
            % Grabs the absolute position of the parent object. 
            % If the parent is a collector node then it is passed on to the
            
            ph = obj.Parent;
            
            %Escape route for when 'get' is called on an empty object
            if isempty(ph)
                val = [0, 0, 0];
                return
            end
            
            %Parent of the collector.
            if isa(ph, 'mvc.model.Collector')
                val = ph.Parent.AbsPosition;
            else
                val = ph.AbsPosition;
            end 
        end
    end
    
    methods % construction / destruction
        
        function obj = ControlSurface(varargin)
            
%             %Check if the user is defining the 'Label'
%             lab = [];
%             ind = find(ismember(varargin, 'Label'), 1);
%             if ~isempty(ind)
%                 lab = varargin{ind + 1};
%                 varargin([ind, ind + 1]) = [];
%             end
            
            %Pass it on           
            obj@awi.model.Component(varargin{:});
            
            %Set the 'Label' value
%             obj.Label = lab;
            
            %Configure collectables
            obj.IsLeafNode = true;
            
            %Extend property groups
            obj.addPropertyGroup('Geometry', ...
                'Eta'   , ['Non-dimensionable spanwise position of the '   , ...
                'start of the control surface with respect to its parent.'], ...
                'xLE'   , ['Normalised chorwise position of the leading '  , ...
                'edge of the control surface with respect to the leading ' , ...
                'edge of the parent.'], ...
                'xTE'  , ['Normalised chorwise position of the trailing '  , ...
                'edge of the control surface with respect to the leading ' , ...
                'edge of the parent.'], ...
                'Max_def' , ['Maximum permitted deflection of the control ', ...
                'surface'], ...
                'Max_rate', ['Maximum permitted rate of deflection of the ', ...
                'control surface']);
            
            %Extend property groups
            obj.addPropertyGroup('Appearance'     , ...
                'FaceColor', 'Colour of the patch', ...
                'FaceAlpha', 'Transparencey of the patch');
            
            %ControlSurface is a primary component in a FE model
            obj.IsFEComponent = true;
            
        end
        
    end
    
    methods % visualisation
        function hg = drawElement(obj, ht, tag)
            
            %TODO - Update the drawing of the control surface so it is in
            %its own reference system. (Allows us to rotate it)
            if nargin < 3
               tag = 'Generic Control Surfaces'; 
            end
            
            %Start with base-class
            %hg = drawElement@awi.model.Component(obj, ht);
            
            %Draw the coordinate system of the hinge
            %   TODO - Update the drawing of the hinge coordinate system
%             hg = drawElement(obj.HingeCoordSys, ht);
            hg = {};
            
            %Position of origin relative to parent
            O = obj.ParentAbsPos;           
            
            %Plot aerodynamic panel outline
            coords = obj.Coords; % use get method once!
            
            %Anything ?
            if isempty(coords)
                
                %NO
                return;
                
            end
            
            %Extend drawing
            xAP = [coords.LE.X ; coords.TE.X] - O(1);
            yAP = [coords.LE.Y ; coords.TE.Y] - O(2);
            zAP = [coords.LE.Z ; coords.TE.Z] - O(3);
            hg{end+1} = surf(xAP, yAP, zAP, 'Parent', ht, ...
                'FaceColor', obj.FaceColor, ...
                'FaceAlpha', obj.FaceAlpha, ...
                'Tag'      , tag); 
        end
    end
    
    methods % converting to FE model
        function FEModel = convertThisToFE(obj, FEModel, StructuralNodes, varargin)
            %convertThisToFE Class specific method for converting an
            %instance of 'awi.model.ControlSurface' into a collection of FE
            %data objects.
            %
            % - Make a grid point at every vertex
            % - Make a coordinate system at the inboard hinge location 
            % - Create a rigid link between the the hinge grid point and 
            %   the other 3 vertices
            % - Make the aerodynamic panels*
            % - Create a aero panel set*
            % - Create a structural set*
            % - Create a spline*
            
            assert(numel(obj) == 1, ['Function only valid for ', ...
                'scalar instance of the ''awi.fe.ControlSurface'' class']);
        
            if nargin < 3
                StructuralNodes = [];
            end
            
            %Start with base-class
            FEModel = convertThisToFE@awi.mixin.FEable(obj, FEModel);
            
            %Grab the vertex coordinates in the global coordinate system
            vertexCoords = obj.Coords;
            coords = [ ...
                vertexCoords.LE.X, vertexCoords.TE.X ; ...
                vertexCoords.LE.Y, vertexCoords.TE.Y ; ...
                vertexCoords.LE.Z, vertexCoords.TE.Z ];
            
            if size(coords, 2) > 4
                return
            end
            
            %Make the aero panels
            AeroPanel = i_generatePanels(obj, vertexCoords);
            
            %Make the nodes associated with these nodes
            g = arrayfun(@(~) awi.fe.Node, 1 : 4, 'Unif', false);
            g = horzcat(g{:});
            set(g, {'X'}, num2cell(coords, 1)');            
            
            %Define the coordinate system
            coordSys = defineHingeCoordSys(obj, vertexCoords);
            hingeCoords = coordSys.A;
            
            %Find the hinge coordinates in the set
            idx = ismember(coords', hingeCoords', 'rows');
            
            %Define the rigid link between the hinge node and the other 3
            Rigid = awi.fe.RigidBar;
            Rigid.NodesI = g(idx);
            Rigid.NodesD = g(~idx)';
            
            %Find closest node along the beam
            r = sqrt(sum((g(idx).X - [StructuralNodes.X]) .^ 2));
            [~, ind] = min(r);
            
            %Add a bush element between the nearest beam node and the
            %control system.
            bush     = awi.fe.BushElement;
            bushprop = awi.fe.BushProp;
            bush.Nodes        = [g(idx) ; StructuralNodes(ind)];
            bush.BushProperty = bushprop;
            
            %Set a high stiffness value for the bush element
            %   - May need to tweak this value depending on the
            %   model. Can give MATRIX/FACTOR DIAGONAL warnings.
            k    = 1e11;
            kNam = arrayfun(@(i) ['K', num2str(i)], 1 : 6, 'Unif', false);
            kVal = repmat({k}, [1, 6]);
            set(bushprop, kNam, kVal);
            
            %Update the FEM
            addFEData(FEModel, g, coordSys, Rigid, AeroPanel, bush, bushprop);
            
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
                switch obj.ParentBeam.SpanVector
                    case 'X'
                        error('Update code for ''SpanVector'' = X');
                    case 'Y'
                        s = abs(diff(vertcat(leCoords.Y), [], 2));
                    case 'Z'
                        s = abs(diff(vertcat(leCoords.Z), [], 2));
                end
                c = max(chrd, [], 2);
                
                %Calculate aero panel size
                [nSpan, nChord] = awi.model.LiftingSurface.getAeroPanelSize(obj, s, c);
                
                %Generate the aerodynamic panels
                AeroPanel = arrayfun(@(~) awi.fe.AeroControlSurf, 1 : numel(cs), 'Unif', false);
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
        function coordSys = defineHingeCoordSys(obj, vertexCoords)
            %defineHingeCoordSys Defines the 'awi.fe.CoordSys' object that
            %sits at the hinge location of the coordinate system.
            
            if nargin < 2
               vertexCoords = obj.Coords; 
            end
            
            %Coordinate where the hinge will be defined?
            opt = {'LE', 'TE'};
            idx = ismember(opt, obj.HingeLine);
            hingeCoords = [ ...
                vertexCoords.(opt{idx}).X(1) ; ...
                vertexCoords.(opt{idx}).Y(1) ; ...
                vertexCoords.(opt{idx}).Z(1)];            
                        
%             %Coordinate that will define the coordinate system axis
%             vectorCoords = [ ...
%                 vertexCoords.(opt{~idx}).X(1) ; ...
%                 vertexCoords.(opt{~idx}).Y(1) ; ...
%                 vertexCoords.(opt{~idx}).Z(1)];
            
            %Need to calculate local z-direction of the aero-plane
            vx = [ ...
                vertexCoords.(opt{~idx}).X(1) - vertexCoords.(opt{idx}).X(1) ; ...
                vertexCoords.(opt{~idx}).Y(1) - vertexCoords.(opt{idx}).Y(1) ; ...
                vertexCoords.(opt{~idx}).Z(1) - vertexCoords.(opt{idx}).Z(1)]; 
            vy = [ ...
                diff(vertexCoords.(opt{idx}).X) ; ...
                diff(vertexCoords.(opt{idx}).Y) ; ...
                diff(vertexCoords.(opt{idx}).Z)];
            vz = cross(vx, vy);
            
            %Now calculate the new x-vector so we can ensure the y-axis is
            %correct - Currently the y-axis will just be aligned with the
            %global axis
            vx = cross(vy, vz);
            
            %Define the coordinate system
            coordSys = awi.fe.CoordSys;
            coordSys.A = hingeCoords;
            coordSys.B = hingeCoords + vz;
%             coordSys.C = vectorCoords;
            coordSys.C = hingeCoords + vx;
            
        end
    end
    
end
