classdef (ConstructOnLoad) Connector < awi.model.Stick
    %Connector Defines a line connection between two components.
    
    properties (AbortSet, SetObservable)
        %Name of the parent object of the tip node (Character Array)
        TipParent
        %Normalised distance along the parent span that the tip node is
        %located at
        TipNormPos
        %Defines the span vector which 'TipNormPos' is defined along.
        TipFlag = 'R';
        %X offset in the global coordinate system of the tip node from the
        %location defined by 'TipNormPos', 'TipParent' & 'TipFlag'
        TipXOffset = 0;
        %Y offset in the global coordinate system of the tip node from the
        %location defined by 'TipNormPos', 'TipParent' & 'TipFlag'
        TipYOffset = 0;
        %Z offset in the global coordinate system of the tip node from the
        %location defined by 'TipNormPos', 'TipParent' & 'TipFlag'
        TipZOffset = 0;
        %Degrees of Freedom (DOFs) that are matched at the root
        RootDOF = '123456';
        %Degrees of Freedom (DOFs) that are matched at the tip
        TipDOF = '123456';
    end    
    
    properties (SetAccess = private)
        %Handle to child connectors
        ChildConnectors
    end
    
    properties %(Dependent)
        %Handle to the parent object of the tip node
        TipParent_
    end
    
    methods % set / get
        function set.TipParent(obj, val)   %set.TipParent  
            %set.TipParent Set method for the property 'TipParent'.
            %
            % 'TipParent' must be a row character array.
            
            validateattributes(val, {'char'}, {'row'}, class(obj), ...
                'TipParent');
            obj.TipParent = val;
        end
        function set.TipNormPos(obj, val)  %set.TipNormPos 
            %set.TipNormPos Set method for the property 'TipNormPos'
            %    - 'TipNormPos' must be a scalar numeric with a value in
            %    the range [0 : 1]
            validateattributes(val, {'numeric'}, {'scalar', ...
                'nonnegative', '<=', 1}, class(obj), 'TipNormPos');
            obj.TipNormPos = val;
        end
        function set.TipFlag(obj, val)     %set.TipFlag    
            %set.TipFlag Set method for the property 'TipFlag'.
            %
            %   - 'TipFlag' must be one of the following tokens:
            %     {'X', 'Y', 'Z', 'R'}.
            
            val = validatestring(val, {'X', 'Y', 'Z', 'R'}, class(obj), ...
                'TipFlag');
            obj.TipFlag =  val;
        end
        function set.TipXOffset(obj, val)  %set.TipXOffset 
            %set.TipXOffset Set method for the property 'TipXOffset'
            %
            %    - 'TipXOffset' must be a scalar numeric
            
            validateattributes(val, {'numeric'}, {'scalar', 'nonnan', ...
                'finite', 'real'}, 'TipXOffset', class(obj));
            obj.TipXOffset = val;
        end
        function set.TipYOffset(obj, val)  %set.TipYOffset 
            %set.TipYOffset Set method for the property 'TipYOffset'
            %
            %    - 'TipYOffset' must be a scalar numeric
            
            validateattributes(val, {'numeric'}, {'scalar', 'nonnan', ...
                'finite', 'real'}, 'TipYOffset', class(obj));
            obj.TipYOffset = val;
        end
        function set.TipZOffset(obj, val)  %set.TipZOffset 
            %set.TipZOffset Set method for the property 'TipZOffset'
            %
            %    - 'TipZOffset' must be a scalar numeric
            
            validateattributes(val, {'numeric'}, {'scalar', 'nonnan', ...
                'finite', 'real'}, 'TipZOffset', class(obj));
            obj.TipZOffset = val;
        end
        function set.TipParent_(obj, val)  %set.TipParent_ 
            %set.TipParent_ Set method for the property 'TipParent_'.
            %
            % 'TipParent_' must be of type, or a subclass of,
            % 'awi.model.Connector'.
            % 'TipParent_' must be a scalar.
            
            if isempty(val)
                obj.TipParent_ = val;
                return
            end
            
            validateattributes(val, {'awi.model.Connector'}, {'scalar'}, ...
                class(obj), 'TipParent_');
            obj.TipParent_ = val;
            
            %Update the list of child connectors in the parent
            addChildConnector(val, obj);
            
        end
        function val = get.TipParent_(obj) %get.TipOffset_
            %get.TipParent_ Get method for the dependent property
            %'TipParent_'.
            
            %This should really be something like
            val = obj.TipParent_;
            if isempty(val)
                %Search the hierachy for the name of the object
                %Get top level ancestor
                ancstr = ancestor(obj, 'awi.model.Entity', 'toplevel');                
                %Search the framework for the object with name 'TipParent'
                val = findall(ancstr, 'Name', obj.TipParent);
                %Assign it so we don't need to traverse the tree again
                if isempty(val)
                    val = [];
                end
                obj.TipParent_ = val;
            end
            
            %Force empty matrix
            if isempty(val) 
                val = [];
            end
            
        end
    end
    
    methods % constructor
        
        function obj = Connector(varargin)
            %Connector Constructor for the class 'awi.model.Connector'.
            
            %Pass it on
            obj@awi.model.Stick(varargin{:});
            
        end
        
    end
    
    methods % helper methods
        function addChildConnector(obj, child)
           %addChildConnector Manages a list of 'awi.model.Connector' 
           %objects that are connected to this object.
            
           %Manage the list...
           if isempty(obj.ChildConnectors)
               %First time through
               obj.ChildConnectors = child;
           elseif ~ismember(obj.ChildConnectors, child)
               %Only add the object if it isn't already in the list
               obj.ChildConnectors = [obj.ChildConnectors, child];  
           end
        end
        function getConnectorXYZ(obj)
            %getConnectorXYZ Calculates the (x,y,z) coordinates of the
            %connector in the parent coordinate system.
            
            %Root position is obvious...
            rootPosition = obj.Position;
            
            %Find parent object of the tip node
            phTip = obj.TipParent_;
            
            %Check if the parent handle has been return correctly
            if isempty(phTip)
                return
            end
            
            %Calculate coordinates of the tip of the stick
            %   - These coordinates are in the coordinate system of the tip
            %   parent object.
            %   - We must convert them to the global system and then into
            %   the coordinate of the root parent object.
            if isa(phTip, 'awi.model.Stick') %Interpolate along the stick
                tipPosition  = s2pos(phTip, obj.TipNormPos, obj.TipFlag);
            else
                tipPosition = phTip.Position;
            end
            tipPosition  = tipPosition + ... %Add offsets in global CS
                [obj.TipXOffset, obj.TipYOffset, obj.TipZOffset];
            tipPosition = tipPosition + ...  %Convert to global CS
                phTip.AbsPosition;
            tipPosition = tipPosition - ...  %Convert to root parent CS
                obj.Parent.AbsPosition;
            
            %Initially define (x,y,z) data as just the root and tip
            obj.XData = [rootPosition(1), tipPosition(1)];
            obj.YData = [rootPosition(2), tipPosition(2)];
            obj.ZData = [rootPosition(3), tipPosition(3)];
        end        
    end
    
    methods % converting to FE
        function FEModel = convertThisToFE(obj, FEModel, varargin)
            %convertThisToFE Converts the 'awi.model.Connector' object to a
            %collection of Finite Element (FE) entities.
            %
            % A connector has the following FE entities:
            %   - Node(s)
            %   - Bush element
            %   - Bush element property
            %   - Joints 
            
            idx = ismember(varargin(1 : 2 : end), 'AddDataToParent');
            if any(idx)
                bAdd2Parent = true;
            else
                bAdd2Parent = false;
            end
            
            %Start at the component level
            FEModel = convertThisToFE@awi.model.Component(obj, FEModel);
            
            %Grab the geometry objects
            RootPar = obj.Parent;
            TipPar  = obj.TipParent_;
            
            if or(isempty(RootPar), isempty(TipPar))
                return
            end
            
            %Find FE model corresponding to the root/tip parent objects
            AllFEM  = flatlist(FEModel.RootModel);
            GObj    = [AllFEM.GeometryObject];            
            RootFEM = AllFEM(ismember(GObj, RootPar));
            TipFEM  = AllFEM(ismember(GObj, TipPar));
            
            %If the FEM is not finished building itself then add some
            %preliminary data and quit out.
            if or(isempty(RootFEM), isempty(TipFEM))
                %Store a reference to the connector and return to this
                %object once the collection has been fully defined.
                FEModel.IsComponent    = true;
                FEModel.GeometryObject = obj;
                %But first we need to add some data so that this FE model
                %is not discarded by the invoking function
                RigidBar = awi.fe.RigidBar;
                addFEData(FEModel, RigidBar);
%                 [Joints, Bush, BushPrp] = awi.fe.FEModel.makeDefaultConnection;
%                 addFEData(FEModel, Joints, Bush, BushPrp);
                return
            end
            
%             %Check that the FE objects that comprise a generic connection
%             %have already been added. If not, add them!
            RigidBar = FEModel.RigidBars;
            if isempty(RigidBar)
                RigidBar = awi.fe.RigidBar;
                addFEData(FEModel, RigidBar);
            end

%             idxJ = arrayfun(@(o) isa(o, 'awi.fe.Joint')      , FEModel.Connections);
%             idxB = arrayfun(@(o) isa(o, 'awi.fe.BushElement'), FEModel.Connections);
%             bNodes   = ~isempty(FEModel.Nodes);
%             bJoint   = any(idxJ);
%             bBush    = any(idxB);
%             bBushPrp = ~isempty(FEModel.BushProps);
%             
%             if bNodes   %Nodes     
%                 Nodes = FEModel.Nodes;
%             else
%                 [Nodes, ~, ~, ~] = awi.fe.FEModel.makeDefaultConnection;
%                 addFEData(FEModel, Nodes);
%             end
%             if bJoint   %Joints    
%                 Joints = FEModel.Connections(idxJ);
%             else
%                 [Joints, ~, ~, ~] = awi.fe.FEModel.makeDefaultConnection;
%                 addFEData(FEModel, Joints);
%             end
%             if bBush    %Bush      
%                 Bush  = FEModel.Connections(idxB);
%             else
%                 [~, ~, Bush, ~] = awi.fe.FEModel.makeDefaultConnection;
%                 addFEData(FEModel, Bush);
%             end
%             if bBushPrp %Bush Prop 
%                 BushPrp = FEModel.BushProps;
%             else
%                 [~, ~, ~, BushPrp] = awi.fe.FEModel.makeDefaultConnection;
%                 addFEData(FEModel, BushPrp);
%             end
            
            %Assign coordinates of root/tip connections to Nodes
            xyzRoot = obj.AbsPosition;
            if isa(TipPar, 'awi.model.Stick')
                xyzTip  = s2pos(TipPar, obj.TipNormPos, obj.TipFlag);
                xyzTip  = xyzTip + TipPar.AbsPosition;
            else
                xyzTip  = TipPar.AbsPosition;
            end
%             set(Nodes, {'X'}, num2cell([xyzRoot ; xyzTip]', 1)');
            
            %Find associated nodes in parent FEM
            RootParNode = i_getConnectorNode(RootFEM, xyzRoot);
            TipParNode  = i_getConnectorNode(TipFEM , xyzTip);
            
            if isempty(RootParNode) || isempty(TipParNode) %Escape route
                return
            end
            
            %Assign nodes to rigid element
            RigidBar.NodesI = RootParNode;
            RigidBar.NodesD = TipParNode;
            RigidBar.CN     = obj.TipDOF;
            
%             %Assign nodes to joints
%             Joints(1).Nodes = [RootParNode ; Nodes(1)];
%             Joints(2).Nodes = [TipParNode  ; Nodes(2)];
%             Joints(1).CB    = obj.RootDOF;
%             Joints(2).CB    = obj.TipDOF;
%             
%             %Assign new nodes to Bush and BushProp to Bush
%             Bush.Nodes        = Nodes';
%             Bush.BushProperty = BushPrp;
            
            %Add FE data to tip parent FEM
            if bAdd2Parent
%                 addFEData(TipFEM, Nodes, Joints, Bush, BushPrp);
                addFEData(TipFEM, RigidBar);
                detach(FEModel);
                delete(FEModel)
                FEModel = [];
            end
                                    
            function Node = i_getConnectorNode(FEModel, xyz)
                %i_getConnectorNode Retrieves the handle to the node in the
                %FEModel 'FEModel' which has coordinates 'xyz'.
                
                Node = [];
                AllNodes = [FEModel.Nodes];
                if isempty(AllNodes)
                    return
                end                
                nodeCoords = [AllNodes.X];
                Node = AllNodes(ismember(nodeCoords', xyz, 'rows'));
            
            end
            
        end        
    end
    
    methods % visualisation
        
        function hg = drawElement(obj, ht, tag)
            %drawElement Class specific method for visualising the
            %'Connector' element.
            
            if nargin < 3
                tag = 'Connectors';
            end
            
            %Set up some defaults
            obj.LineWidth       = 1;
            obj.LineColor       = 'r';
            obj.MarkerStyle     = 'o';
            obj.MarkerFaceColor = 'r';
            obj.MarkerEdgeColor = 'k';
            
            %Define the (x,y,z) data - In parent coordinate system
            getConnectorXYZ(obj);
            
            %Need to define 'XData', 'YData', 'ZData' in local coord system
            pos = obj.Position;
            obj.XData = obj.XData - pos(1);
            obj.YData = obj.YData - pos(2);
            obj.ZData = obj.ZData - pos(3);

            %Pass it on
            hg = drawElement@awi.model.Stick(obj, ht, tag);
            
        end
        
    end
    
end

