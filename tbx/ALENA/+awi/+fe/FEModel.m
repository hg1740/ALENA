classdef FEModel < awi.fe.FEBaseClass & matlab.mixin.CustomDisplay
    %FEModel Defines a collection of standard Finite Element entities.
    %
    % TODO - Change all FEProps to dynamic properties for a more
    % lightweight implementation. Q - How do heterogeneous objects work in
    % this instance? Only display commen properties?
    % TODO - Add a custom display so we only display properties that are
    % non-empty.
    % TODO - Add a method for adding properties/bulk data types in the
    % constructor. (Think MVC proeprty sets)
    % TODO - Add property which flags whether this model is a residual
    % structure (i.e. massid = ... or BEGIN SUPER...)
    %
    
    %Primary FE entities that can be part of a FEModel
    properties (SetAccess = private)
        %Coordinate Systems
        CoordSys
        %Nodes
        Nodes
        %Scalar Points
        ScalarNodes
        %Point Masses
        PointMasses
        %Beam Elements
        Beams
        %Beam Element Properties
        BeamProps
        %Bush Element Properties
        BushProps
        %Scalar Mass Elements
        ScalarMasses
        %Scalar Mass Element Properties
        ScalarMassProps
        %Damper Elements
        Dampers
        %Damper Element Properties
        DamperProps
        %Spring Elements
        Springs
        %Spring Elemet Properties
        SpringProps
        %Materials
        Materials
        %Rigid Bars/Connections
        RigidBars
        %Aerodynamic Panels
        AeroPanels
        %Aerodynamic Control Surfaces
        AeroControlSurf
        %Aerodnamic Panel Properties
        AeroProps
        %Structural Sets
        StructuralSets
        %Aerodynamic Sets
        AeroSets
        %Aeroelastic Splines
        Splines
        %Connections between the FE model and its parent
        Connections
        %Mass Sets
        MassGroups
        %Loads applied at a point in the model
        PointLoads
        %DMIs
        DMI
    end
    
    %Optimisation entities that can be part of a FEModel
    properties (SetAccess = private)
        %Design variables
        DesVar = awi.fe.opt.DesignVariable.empty;
        %Design property relations
        DesPropRel = awi.fe.opt.PropRelation.empty;
        %Design propery relations (using equations)
        DesPropRelEqn = awi.fe.opt.PropRelationEqn.empty;
        %Design responses
        DesResp = awi.fe.opt.DesignResponse.empty;
        %Design response equations
        DesRespEqn = awi.fe.opt.DesignResponseEqn.empty;
        %Design constraints
        DesConstr = awi.fe.opt.DesignConstraint.empty;
        %Sets of design constraints
        DesConstrSet = awi.fe.opt.DesignConstraintSet.empty;
        %Design constants
        DesTable = awi.fe.opt.DesignTable.empty;
        %Design equations
        DesEqn = awi.fe.opt.DesignEquation.empty;
    end
    
    %Collections of 'awi.fe' objects that can be part of a FEModel
    properties (SetAccess = private)
        NodeCollection   
        PanelCollection
        PanelPropCollection
    end
    
    %Identifiers
    properties
        %Name of this FEModel
        Name = '';
        %Determines whether the FEModel is a component in an assembly of
        %other FEModels.
        IsComponent = false;
        %Parent of the FE model
        Parent
        %Children of the FE model
        Children
        %Handle to the geoemtry object that this component is derived from
        GeometryObject            
    end
    
    %Helper properties
    properties
        %Function for logging the progress of the FE conversion process
        LogFcn = @(s) fprintf('%s\n', s);
        %Level in the hierachy (used during logging of FE conversion)
        HierarchyLevel = 0;
    end
    
    %Helper properties
    properties (SetAccess = private)
        %Maximum ID number in this collection
        MaxIDNumber
        %List of parts in the model
        PartList
    end
    
    %Helper properties
    properties (Dependent, Hidden = true)
        %FE model at the top of the hierachy
        RootModel
        %Nodes belonging to beam elements
        BeamNodes
        %Indentation character
        IndentChar
    end
    
    %Helper properties
    properties (Constant, Hidden = true)
        %Character for denoting a new level in the hierachy
        IndentTok = '  ';
        %Mapping between FE objects and FEModel properties        
        map = { ...
            'awi.fe.CoordSys'         , 'CoordSys'       , 'CORD2R'  ; ...
            'awi.fe.Node'             , 'Nodes'          , 'GRID'    ; ...
            'awi.fe.ScalarPoint'      , 'ScalarNodes'    , 'SPOINT'  ; ...
            'awi.fe.Beam'             , 'Beams'          , 'CBEAM'   ; ...
            'awi.fe.BeamProp'         , 'BeamProps'      , 'PBEAM'   ; ...
            'awi.fe.BeamCrossSection' , 'BeamProps'      , 'PBEAML'  ; ...
            'awi.fe.Material'         , 'Materials'      , 'MAT1'    ; ...
            'awi.fe.RigidBar'         , 'RigidBars'      , 'RBE2'    ; ...
            'awi.fe.AeroPanel'        , 'AeroPanels'     , 'CAERO1'  ; ...
            'awi.fe.AeroProp'         , 'AeroProps'      , 'PAERO'   ; ...
            'awi.fe.AeroControlSurf'  , 'AeroControlSurf', 'AESURF'  ; ...
            'awi.fe.StructuralSet'    , 'StructuralSets' , 'SET1'    ; ...
            'awi.fe.AeroPanelSet'     , 'AeroSets'       , 'AELIST'  ; ...
            'awi.fe.AeroelasticSpline', 'Splines'        , 'SPLINE7' ; ...
            'awi.fe.PointMass'        , 'PointMasses'    , 'CONM2'   ; ...
            'awi.fe.Joint'            , 'Connections'    , 'RJOINT'  ; ...
            'awi.fe.BushElement'      , 'Connections'    , 'CBUSH'   ; ...
            'awi.fe.BushProp'         , 'BushProps'      , 'PBUSH'   ; ...
            'awi.fe.ScalarMass'       , 'ScalarMasses'   , 'CMASS'   ; ...
            'awi.fe.ScalarMassProp'   , 'ScalarMassProps', 'PMASS'   ; ...
            'awi.fe.Damper'           , 'Dampers'        , 'CDAMP1'  ; ...
            'awi.fe.DamperProp'       , 'DamperProps'    , 'PDAMP'   ; ...
            'awi.fe.Spring'           , 'Springs'        , 'CELAS'   ; ...
            'awi.fe.SpringProp'       , 'SpringProps'    , 'PELAS'   ; ...
            'awi.fe.MassSet'          , 'MassGroups'     , 'MASSSET' ; ...
            'awi.fe.PointLoad'        , 'PointLoads'     , 'FORCE'   ; ...  %ACTUALLY IT IS FORCE/MOMENT!!
            'awi.fe.opt.DesignVariable'     , 'DesVar'       , 'DESVAR'  ; ...
            'awi.fe.opt.PropRelation'       , 'DesPropRel'   , 'DVPREL1' ; ...
            'awi.fe.opt.PropRelationEqn'    , 'DesPropRelEqn', 'DVPREL2' ; ...
            'awi.fe.opt.DesignResponse'     , 'DesResp'      , 'DRESP1'  ; ...
            'awi.fe.opt.DesignResponseEqn'  , 'DesRespEqn'   , 'DRESP2'  ; ...
            'awi.fe.opt.DesignConstraint'   , 'DesConstr'    , 'DCONSTR' ; ...
            'awi.fe.opt.DesignConstraintSet', 'DesConstrSet' , 'DCONADD' ; ...
            'awi.fe.opt.DesignTable'        , 'DesTable'     , 'DTABLE'  ; ...
            'awi.fe.opt.DesignEquation'     , 'DesEqn'       , 'DEQATN'  ; ...
            'awi.fe.NodeCollection'         , 'NodeCollection', 'GRID'   ; ...
            'awi.fe.PanelCollection'        , 'PanelCollection', 'CQUAD4'; ...
            'awi.fe.PanelPropCollection'    , 'PanelPropCollection', 'PSHELL'; ...
            'awi.fe.DMI'                    , 'DMI'         ,'DMI'};
        
    end
    
    methods % set / get
        function set.Name(obj, val)           %set.Name
            %set.Name Set method for the property 'Name'.
            %
            % 'Name' must be a row vector of characters.
            validateattributes(val, {'char'}, {'row'}, class(obj), 'Name');
            obj.Name = val;
        end
        function set.IsComponent(obj, val)    %set.IsComponent
            %set.IsComponent Set method for the property 'IsComponent'.
            %
            % 'IsComponent' must be a scalar logical.
            validateattributes(val, {'logical'}, {'scalar'}, ...
                class(obj), 'IsComponent');
            obj.IsComponent = val;
        end
        function set.Parent(obj, val)         %set.Parent
            %set.Parent Set method for the property 'Parent'.
            %
            % 'Parent' must be a scalar instane of the 'awi.fe.FEModel'
            % class.
            
            if isempty(val) %Always okay
                obj.Parent = [];
                return
            end
            
            validateattributes(val, {'awi.fe.FEModel'}, {'scalar'}, ...
                class(obj), 'Parent');
            obj.Parent = val;
        end
        function set.Children(obj, val)       %set.Children
            %set.Children Set method for the property 'Children'.
            %
            % 'Children' must be a row vector of the 'awi.fe.FEModel'
            % class.
            
            if isempty(val) %Always okay
                obj.Children = [];
                return
            end
            
            validateattributes(val, {'awi.fe.FEModel'}, {'row'}, ...
                class(obj), 'Children');
            obj.Children = val;
        end
        function set.GeometryObject(obj, val) %set.GeometryObject
            %set.GeometryObject Set method for the property
            %'GeometryObject'.
            %
            % 'GeometryObject' must be a scalar object which derives from
            % the 'awi.model.Component' object.
            
            validateattributes(val, {'awi.model.Component'}, {'scalar'}, ...
                class(obj), 'GeometryObject');
            obj.GeometryObject = val;
        end
        function set.LogFcn(obj, val)         %set.LogFcn
            assert(isa(val, 'function_handle'), ['Expected ''LogFcn'' ', ...
                'to be a function handle.']);
            obj.LogFcn = val;
        end
        function set.HierarchyLevel(obj, val) %set.HierarchyLevel
           validateattributes(val, {'numeric'}, {'scalar', 'integer', 'nonnegative'}, ...
               class(obj), 'HierarchyLevel');
           obj.HierarchyLevel = val;
        end
        function val = get.RootModel(obj)     %get.RootModel
            %get.RootModel Get method for the dependent property
            %'RootModel'.
            %
            % 'RootModel' is the 'awi.fe.FEModel' at the top of the
            % hierachy.
            
            %Sensible start
            val = obj;
            
            %Simply go up the hierachy until we are at the top
            while ~isempty(val.Parent)
                val = val.Parent;
            end
            
        end
        function val = get.BeamNodes(obj)     %get.BeamNodes
            %get.BeamNodes Get method for the dependent property
            %'BeamNodes'.
            
            val = [];
            if isempty(obj.Nodes) || isempty(obj.Beams) %Escape route
                return
            end
            val = [obj.Beams.Nodes];
            
        end
        function val = get.IndentChar(obj)    %get.IndentChar
            val = repmat(obj.IndentTok, [1, obj.HierarchyLevel]);
        end
    end
    
    methods % construction
        function obj = FEModel
            %FEModel Constructor for the 'awi.fe.FEModel' class.
            %
            % Actions performed:
            %   - Make a note of the names of the finite element properties
            %     of this object.
            %
            
            %Make a note of the property names
            addFEProp(obj       , 'Children'       , 'CoordSys'     , 'Nodes'          , 'ScalarNodes', ...
                'PointMasses'   , 'Beams'          , 'BeamProps'    , 'BushProps'      , ...
                'ScalarMasses'  , 'ScalarMassProps', 'Springs'      , 'SpringProps'    , ...
                'Dampers'       , 'DamperProps'    , 'Materials'    , 'RigidBars'      , ...
                'AeroPanels'    , 'AeroControlSurf', 'AeroProps'    , 'StructuralSets' , ...
                'AeroSets'      , 'Splines'        , 'Connections'  , 'MassGroups'     , ...
                'PointLoads'    , 'DesVar'         , 'DesPropRel'   , 'DesPropRelEqn'  , ...
                'DesResp'       , 'DesRespEqn'     , 'DesConstr'    , 'DesConstrSet'   , ...
                'DesTable'      , 'DesEqn'         , 'DMI', ...
                'NodeCollection', 'PanelCollection', 'PanelPropCollection');
            
        end
    end
    
    methods (Sealed) % handling a collection of FE models
        function addFEData(obj, varargin)
            %addFEData Adds a finite element entity to the model.
            
            %Check the input is valid
            idx_ = cellfun(@(x) isa(x, 'awi.fe.FEBaseClass'), varargin);
            
            %Only retain valid FE entities
            varargin = varargin(idx_);
            if isempty(varargin)
                return
            end
            
            %What have we got?
            allClass    = cellfun(@class, varargin, 'Unif', false);
            uniqueClass = unique(allClass);
            
            %Catch any heterogeneous inputs
            idxHetero = contains(uniqueClass, 'awi.fe.FEBaseClass');
            if any(idxHetero)
                
                %Strip from the main list
                HeteroFE = varargin(idxHetero);
                varargin(idxHetero) = [];
                
                %Determine unique subclasses of 'awi.fe.FEBaseClass'.
                HeteroFE    = num2cell(horzcat(HeteroFE{:}));
                heteroClass = cellfun(@class, HeteroFE, 'Unif', false);
                uCls        = unique(heteroClass);
                
                %Add the heterogeneous data
                i_addFEData(obj, uCls, heteroClass, HeteroFE);                
     
            end
            
            if isempty(varargin) %Escape route
                return
            end
            
            %Add the non-heterogeneous data
            i_addFEData(obj, uniqueClass, allClass, varargin);
            
            function i_addFEData(obj, uniqueClass, allClass, FEData)
                
                %Add each
                for i = 1 : numel(uniqueClass)
                    
                    %Down-select the type by interrogating the object class
                    cls = uniqueClass{i};
                    
                    %Find these objects in the map and determine the associated
                    %property name
                    idx     = ismember(obj.map(:, 1), cls);
                    if ~any(idx) %Escape route
                        continue
                    end
                    prpName = obj.map{idx, 2};
                    
                    %Grab the items from 'varargin' & force row vector
                    idx   = ismember(allClass, cls);
                    items = cellfun(@(x) x(:)', FEData(idx), 'Unif', false);
                    items = horzcat(items{:});
                    
                    %Add or append the item to the list
                    % - Let the set method handle the details...
                    if isempty(obj.(prpName))                        
                        obj.(prpName) = items;
                    else
                        obj.(prpName) = [obj.(prpName), items];
                    end
                    
                end
                
            end
            
        end
        function addPart(obj, nam, FEData)
            %addPart Assigns a collection of 'awi.fe' data to a part and
            %stores a reference of it.
            
            if isempty(nam) || isempty(FEData)
                return
            end
            
            assert(isvarname(nam), 'The part name must be a valid variable name');
            validateattributes(FEData, {'awi.fe.FEBaseClass'}, {'nonempty'}, ...
                class(obj), 'FEData');
            
            if isempty(obj.PartList)
                obj.PartList = {nam, FEData};
            else
            	obj.PartList = [obj.PartList ; {nam, FEData}];
            end
            
        end
        function lst = flatlist(obj)
            %flatlist Returns a list of all the 'awi.fe.FEModel' objects
            %and their children.
            %
            % Similar to the 'flatlist' method of 'mvc.mixin.Collectable'.
            
            if isempty(horzcat(obj.Children))
                
                %No children so just return the object
                lst = obj;
                
            else
                
                %Return this object and any children in a flat list
                lst = horzcat(obj, flatlist(horzcat(obj.Children)));
                
            end
            
        end
        function varargout = detach(obj)
            %detach Disconnects obj from its parent gracefully, by first
            %removing obj from list of children in parent, and then
            %removing obj's parent
            
            %For each supplied object
            for i = 1:numel(obj)
                
                %Get parent
                par = obj(i).Parent;
                
                %Anything ?
                if isempty(par) || isempty(par.Children)
                    continue;
                end
                
                %Find this object in parent's children
                [b,idx] = ismember(obj(i), par.Children);
                
                %Anything ?
                if ~b
                    warning('unexpected connectivity error'); % or maybe error ?
                    continue;
                end
                
                %Break link
                par.Children(idx) = [];
                obj(i).Parent = [];
                
            end
            
            if nargout == 1
                varargout{1} = obj;
            end
            
        end
        function assignIDnumbers(obj, ExtraFEObj)
            %assignIDnumbers Assignes a unique ID number to every 'awi.fe'
            %object in the 'FEModel'.
            
            if isempty(obj) %Escape route
                return
            end
            
            if nargin < 2
                ExtraFEObj = [];
            else
                validateattributes(ExtraFEObj, {'awi.fe.FEBaseClass'}, ...
                    {'row'}, class(obj), 'ExtraFEObj');
            end
            
            %Grab names of the FE-objects in the model
            pNames = obj(1).map(:, 2);
            
            %Treat 'awi.fe.AeroPanel' differently
            idxAP = ismember(obj(1).map(:, 1), 'awi.fe.AeroPanel');
            pNames(idxAP) = [];
            
            %Grab FE-objects
            FE_Obj = get(obj, pNames);
            
            %Stack the data into a heterogeneous array
            FE_Obj = FE_Obj(:);
            FE_Obj = horzcat(FE_Obj{:}, ExtraFEObj);
            
            %Generate and assign ID numbers - TODO: Take minimum of ID0
            %numbers
            nObj = numel(FE_Obj);
            id   = obj(1).ID0 : obj(1).ID0 + (nObj - 1);
            set(FE_Obj, {'ID'}, num2cell(id)');
            [obj.MaxIDNumber] = deal(id(end));
            
            %Assign ID numbers for any 'awi.fe.MassModel' objects
            MassModel = obj(arrayfun(@(o) isa(o, 'awi.fe.MassModel'), obj));
            if ~isempty(MassModel)
                mid0 = MassModel(1).MassID0;
                mid  = mid0 :  mid0 + numel(MassModel) - 1;
                set(MassModel, {'ID'}, num2cell(mid)');
            end
            
            %Grab the 'awi.fe.AeroPanel' objects
            AP = get(obj, obj(1).map(idxAP, 2));
            AP(cellfun(@isempty, AP)) = [];
            if isempty(AP)
                return
            end
            AP = horzcat(AP{:});
            
            %If there are no structural 'awi.fe' objects defined then we
            %can get a situation where the current max ID is 0
            if isempty(id)
                id = AP(1).ID0;
            end
            
            %Assign the 'awi.fe.AeroPanel' ID numbers            
            nP  = cumsum([AP.NumPanels]);
            ida = (id(end) + 1) + nP;
            ida = [id(end) + 1, ida(1 : end - 1)];
            set(AP, {'ID'}, num2cell(ida)');
            
            %Update the current maximum ID number
            %   - TODO: This should be "set(obj, 'MaxIDNumber', ida(end));"
            %   but because MaxIDNumebr is 'SetAccess = private' we can't
            %   set it in the superclass 'set' method. Seems a bit silly...
            for i = 1 : numel(obj)
                obj(i).MaxIDNumber = ida(end);
            end
            
            %TODO - Assign the row/col numbers for the W2GJ entries
            
        end
        function updateConnectors(obj)
            %updateConnectors Updates any FE models that are derived from
            %'awi.model.Connector' objects.
            
            FEModels = flatlist(obj);
            
            %Look for 'Connector' derived FE models
            idxC = arrayfun(@(fem) strcmp(class(fem.GeometryObject), ...
                'awi.model.Connector'), FEModels); %#ok<STISA>
            if ~any(idxC)%Escape route
                return
            end
            
            %Update the FE data in these models
            ConModels = FEModels(idxC);
            ConModels = num2cell(ConModels);
            for ii = 1 : numel(ConModels)
                ConModels{ii} = convertThisToFE( ...
                    ConModels{ii}.GeometryObject, ConModels{ii}, ...
                    'AddDataToParent', true);
            end
            
        end
        function connectParentAndChild(obj)
            %connectParentAndChild Creates a connection between the FE
            %model (the parent) and any child FE models.
            
            %For an array of FE models just do each in turn
            if numel(obj) > 1
                arrayfun(@(i) connectParentAndChild(obj(i)), ...
                    1 : numel(obj), 'Unif', false);
                return
            end
            
            %For now, just make sure we only deal with one parent FE model
            assert(numel(obj) == 1, 'Expected the number of parent FE models to be equal to 1.');
            
            %No point continuing if the FE model doesn't have any
            %children...
            if isempty(obj.Children)
                return
            end
            
            %Grab the parent nodes
            parNodes = i_getBeamNodes(obj);
            if isempty(parNodes)
                parCoords = [];
            else
                %Grab coordinates of any nodes belonging to the parent
                parCoords = [parNodes.X];
                xPar = parCoords(1, :);
                yPar = parCoords(2, :);
                zPar = parCoords(3, :);
            end
            
            
            
            %For each child in turn...
            for iCh = 1 : numel(obj.Children)
                
                %Which child?
                ch = obj.Children(iCh);
                
                %Grab the nodes
                chNodes = i_getBeamNodes(ch);
                if isempty(chNodes) || ~isempty(ch.Connections) || isempty(parCoords) %Escape route
                    %Pass it on
                    connectParentAndChild(obj.Children(iCh));
                    continue
                end
                
                %Grab coordinates
                chCoords  = [chNodes.X];
                xCh = chCoords(1, :)';
                yCh = chCoords(2, :)';
                zCh = chCoords(3, :)';
                
                %Is the underlying geometry object defined using the
                %'Hierachy Set' ('hSet')? If so, treat differently...
                if strcmp(ch.GeometryObject.ActiveSet, 'hSet') 
                    
                    %When using the 'hSet' we want to force the root of the
                    %child to join to the root of the parent
                    i_generateConnections(ch, parNodes, chNodes(1), ...
                        xPar, yPar, zPar, xCh(1), yCh(1), zCh(1), ...
                        ch.GeometryObject.RootDOF);
                    
                    %Can't make any more connections of the tip parent has
                    %not been defined...
                    if ~isempty(ch.GeometryObject.TipParent)
                        
                        %Find the FEM corresponding to the tip parent
                        tipPar = ch.GeometryObject.TipParent_;
                        models = flatlist(obj.RootModel);
                        idx    = ismember([models.GeometryObject], tipPar);
                        tipFEM = models(idx);
                        
                        if ~isempty(tipFEM)
                            
                            %Get the ndoes
                            tipNodes = i_getBeamNodes(tipFEM);
                            
                            %Grab the coordinates
                            coords = [tipNodes.X];
                            xTip = coords(1, :);
                            yTip = coords(2, :);
                            zTip = coords(3, :);
                            
                            %Make the connections
                            i_generateConnections(ch, tipNodes, chNodes(end), ...
                                xTip, yTip, zTip, xCh(end), yCh(end), ...
                                zCh(end), ch.GeometryObject.TipDOF);
                            
                        end
                        
                    end
                    
                else
                    
                    %Make the connections - Assume fully-fixed
                    i_generateConnections(ch, parNodes, chNodes, xPar, ...
                        yPar, zPar, xCh, yCh, zCh);
                    
                end
                
                %Pass it on
                connectParentAndChild(obj.Children(iCh));
                
            end
            
            function nodes = i_getBeamNodes(FEM)
                %getBeamNodes Retrieves the nodes belonging to all the beam
                %elements in the model. If no beams are defined it just
                %returns the nodes in the model.
                
                %We actually want to use the nodes belonging to the beam
                %elements as we can gurantee that these DOFs will remain in
                %the 'a_set' matrices of Nastran. Also, they actually
                %belong to the structure, as opposed to the planform nodes
                %which are attached via rigid connections.
                %   - Unless of course the FE model doesn't have any beam
                %   elements in which case just go with the nodes...
                if isempty(FEM.Beams)
                    nodes = FEM.Nodes;
                else
                    %Grab the GA/GB nodes for the beams
                    nodes = [FEM.Beams.Nodes];
                    if isempty(nodes)
                        return
                    end
                    %The beam elements should have been defined along the
                    %span of the beam so we can assume there are duplicate
                    %nodes in the object array - Be selective and take all
                    %GA and only the last GB node.
                    nodes = [nodes(1, :), nodes(2, end)];
                end
                
            end
            
            function i_generateConnections(ch, parNodes, chNodes, xPar, yPar, zPar, xCh, yCh, zCh, dof)
                %i_generateConnections
                
                if nargin < 10 %Assume fully-fixed connection
                    dof = 123456;
                end
                
                %Determine the minimum distance
                r = sqrt((xCh - xPar).^2 + (yCh - yPar).^2 + (zCh - zPar).^2);
                
                %Assume that any node less than 'tol' away from a parent
                %node must be coincident
                tol = 1e-6;
                bCoincident = r < tol;
                
                %Can only handle the case where we have a single coincident
                %node - TODO: Update this.
                if nnz(bCoincident) > 1
                    error(['Unable to handle the case where a parent ', ...
                        'FE model has more than one coincident node ', ...
                        'with one of its children']);
                end
                
                %Index into the parent and child dimensions
                if numel(chNodes) > 1
                    idxPar = any(bCoincident);
                else
                    idxPar = bCoincident;
                end
                idxCh  = any(bCoincident');
                
                if any(idxCh) %Any coincident?
                    %If we have coincident nodes then create a joint
                    %between them...
                    
                    %Which nodes are we connecting?
                    n = [parNodes(idxPar) ; chNodes(idxCh)];
                    
                    %Make the joint
                    j = awi.fe.Joint;
                    
                    %Always assume the constrained DOFs are in the child
                    %coordinate system!
                    j.Nodes = n;
                    j.CB    = dof;
                    
                    %If a FWT is used then we need to adjust the output
                    %coordinate system of the child node to yield the
                    %behaviour described by the flare angle
                    if isa(ch.GeometryObject, 'awi.model.LiftingSurface') && ...
                            ch.GeometryObject.IsFoldingWingTip
                        
                        switch ch.GeometryObject.HingeType
                            case 'pinned'
                                j.CB = 12346;
                            case 'fully-fixed'
                                j.CB = 1234546;
                            case 'ball-socket'
                                j.CB = 123;
                        end
                        
                        %Make a new coordinate system for the flare angle
                        [~, flareRotMatrix] = grabHingeCoordSys(ch.GeometryObject);
                        
                        cs   = awi.fe.CoordSys;
                        cs.A = n(2).X;
                        cs.B = cs.A + flareRotMatrix(:, 3);
                        cs.C = cs.A + flareRotMatrix(:, 1);
                        n(2).OutputCoordSys = cs;
                        
                        %Add bushing element with custom stiffness
                        [~, ~, b, bp] = awi.fe.FEModel.makeDefaultConnection;
                        b.Nodes = n;
                        b.CoordSys = cs;
                        b.BushProperty = bp;                        
                        pv = [{'K1', 'K2', 'K3', 'K4', 'K5', 'K6'} ;
                            num2cell(ch.GeometryObject.HingeStiffness)];                        
                        set(bp, pv{:});
                        
                        addFEData(ch, cs, b, bp);              
                    end
                    
                    %Add the joint to the child FE model
                    addFEData(ch, j);
                    
                else
                    %If we have no coincident nodes then create a
                    %connection between the two closest nodes...
                    
                    %Find the minimum distance and the parent/child nodes
                    %this corresponds to
                    rMin = min(min(r));
                    rInd = find(r == rMin, 1);
                    [chInd, parInd] = ind2sub(size(r), rInd);
                    pNode = parNodes(parInd);
                    cNode = chNodes(chInd);
                    
                    %                     %Connect the parent/child nodes
                    %                     if dof == 123456
                    %                         %If the connection is fully fixed then connect using an
                    %                         %RBE2 element.
                    %                         r = awi.fe.RigidBar;
                    %                         r.NodesI = pNode;
                    %                         r.NodesD = cNode;
                    %                         r.CN     = 123456;
                    %                     else
                    %                         %If a DOF has been released then create a new node
                    %                         %at the child node location and connect the node
                    %                         %using a joint.
                    %
                    %                     end
                    
                    %For each node create a node and a joint, then add a
                    %bushing element between them.
                    %   n(1) = Coincident node at the parent
                    %   n(2) = Coincident node at the child
                    [n, j, b, bp] = awi.fe.FEModel.makeDefaultConnection;
                    
                    %Assign data
                    %   - Same (x,y,z) position
                    %   - Same coordinate system (enables correct
                    %   orientation of the joint)
                    set(n, {'X'}, get([pNode, cNode], 'X'));
                    set(n, {'OutputCoordSys'}, get([pNode, cNode], 'OutputCoordSys'));
                    set(j(1), 'Nodes', [pNode ; n(1)], 'CB', 123456);
                    set(j(2), 'Nodes', [cNode ; n(2)], 'CB', dof);
                    b.Nodes = n';
                    b.BushProperty = bp;
                    
                    %Add it to the child FE model
                    addFEData(ch, n, j, b, bp);
                    
                end
                
            end
            
        end
    end
    
    methods % exporting/writing FE data to a file
        function varargout = export(obj, dn, varargin)
            %export Exports the FE models in MSC.Nastran file format to a
            %user specified folder.
            %
            % Parameter inputs:
            %
            %   * 'CollectChildModels' : If the value of this parameter is
            %   set to 'true' then all of the data for ALL of the FE models
            %   in the collection will be written into a single file.
            %   Default = false
            %
            %   * 'AddHelpfulComments' : If the value of this parameter is
            %   set to 'true' then comments will be provided in the bulk
            %   data files which describe what the data is and what the
            %   numbering scheme is for the different FE components.
            %   Default = true
            %
            %   * 'SplitByBulkType' : If the value of this parameter is set
            %   to 'true' then the data in the model is split by the type
            %   of bulk data entry instead of by parent/child.
            %   Default = false
            %
            %   * 'FullyQualifiedPath' : If the value of this parameter is
            %   set to 'true' then the filepaths in the header file will be
            %   'fully-qualified', meaning that they will be the full
            %   file path and not just relative to the folder in which they
            %   set. (i.e. a local file path)
            %   Default = false
            %
            %   * 'WriteHeaderFile' : If the value of this parameter is set
            %   to 'true' then a header file will be written which
            %   includes the 'Executive Control', 'Case Control' and 'Bulk
            %   Data' sections of the MSC.Nastran data file.
            %
            %   * 'WriteOptimisationData' : If the value of this parameter
            %   is set to 'true' then any 'awi.fe.opt' objects in the
            %   collection will be written as bulk data entries in the
            %   '.bdf' files.
            %   Default = true
            %
            %   * 'WriteDesignModel' : If the value of this parameter is
            %   set to 'true' then any 'awi.fe.opt.DesVar' objects and
            %   associated 'aw.fe.DesignPropertyRelation' objects are
            %   written into a seperate file title
            %   "<Name>_DesignModel.bdf".
            %   Default = false
            
            varargout = {{''}, {''}};
            
            %Cannot proceed if no data has been specified
            bData = arrayfun(@(fe) fe.HasFEData, obj);
            if ~any(bData)
                return
            end
            
            if nargin < 2 || isempty(dn) %Directory provided?
                %Ask the user
                dn = uigetdir(pwd, ['Select a folder for the FEM ', ...
                    'bulk data files to be written in']);
                %Valid?
                if isnumeric(dn)
                    return
                end
            else
                assert(isfolder(dn), ['The export function must be ', ...
                    'given a valid directory to write the FEM '  , ...
                    'bulk data files.']);
            end
            
            %Parse inputs
            p = inputParser;
            val_func = @(x)validateattributes(x, {'logical'}, {'scalar'});
            addRequired(p, 'dn', @isfolder);
            addParameter(p, 'CollectChildModels'   , false, val_func);
            addParameter(p, 'AddHelpfulComments'   , true , @(x)validateattributes(x, {'logical'}, {'scalar'}));
            addParameter(p, 'SplitByBulkType'      , false, @(x)validateattributes(x, {'logical'}, {'scalar'}));
            addParameter(p, 'FullyQualifiedPath'   , false, @(x)validateattributes(x, {'logical'}, {'scalar'}));
            addParameter(p, 'Solution'             , 103  , @(x)validateattributes(x, {'numeric'}, {'scalar', 'integer'}));
            addParameter(p, 'NModes'      , 20, @(x)validateattributes(x, {'numeric'}, {'integer', 'scalar'} , 'writeNastranHeaderFile', 'NModes'));
            addParameter(p, 'MaxFreq'     , [], @(x)validateattributes(x, {'numeric'}, {'scalar', 'positive'}, 'writeNastranHeaderFile', 'MaxFreq'));
            addParameter(p, 'WriteHeaderFile'      , true , @(x)validateattributes(x, {'logical'}, {'scalar'}));
            addParameter(p, 'WriteOptimisationData', true , @(x)validateattributes(x, {'logical'}, {'scalar'}));
            addParameter(p, 'WriteDesignModel'     , false, @(x)validateattributes(x, {'logical'}, {'scalar'}));
            parse(p, dn, varargin{:});
            
            %Collect all the models into a single file?
            if p.Results.CollectChildModels
                fem = flatlist(obj);
            else
                fem = obj;
                % TODO ...
                %Write the .xml file describing the model numbering scheme
                % ModelNumberingScheme.xml // ModelNumberingScheme.h5
            end
            
            %Write the data to a file
            partNames = writeToFile(fem, dn    , ...
                p.Results.AddHelpfulComments   , ...
                p.Results.CollectChildModels   , ...
                p.Results.WriteOptimisationData, ...
                p.Results.WriteDesignModel);
            
            if isempty(partNames) %Escape route
                return
            end
            
            %Fully-qualified file names?
            if p.Results.FullyQualifiedPath
                partNames = fullfile(dn, partNames);
            end
            
            %Write the header file
            if p.Results.WriteHeaderFile
                hName = fullfile(dn, 'NastranHeaderFile.dat');
                awi.fe.FEBaseClass.writeNastranHeaderFile(hName, ...
                    'Solution'    , p.Results.Solution, ...
                    'IncludeFiles', partNames, ...
                    'NModes'      , p.Results.NModes, ...
                    'MaxFreq'     , p.Results.MaxFreq);
            else
                hName = '';
            end
            
            %Define outputs
            if nargout > 0
                varargout{1} = hName;
            end
            if nargout > 1
                varargout{2} = partNames;
            end
            
        end
        function fNam = writeToFile(obj, dn, bComment, bCollect, bWriteOpt, bDesModel, bMassCase)
            %writeToFile Writes the model to a file in the folder specified
            %by 'dn'.
            %
            % Inputs:
            %
            %   * 'dn' : Directory location where the FE model bulk data
            %   files will be written.
            %   * 'bComment'  : Defines whether the model summary comments
            %     are written into the file.
            %   * 'bCollect'  : Defines whether the collection of FE models
            %     will be written into a single file or whether each FE
            %     model in the hierarchy will be written into its own file.
            %   * 'bWriteOpt' : Defines whether any optimisation data is
            %     written into the file.
            %
            % Optional inputs:
            %
            %   * 'bDesModel' : Defines whether any optimisation data
            %     and/or designed properties are to be written into a
            %     seperate file.
            %   * 'bMassCase' : Defines whether this FE model represents a
            %     new 'Mass Case', or in MSC.Nastran parlance a "residual
            %     structure".
            
            if nargin < 6 || isempty(bDesModel)
                bDesModel = false;
            end
            if nargin < 7 || isempty(bMassCase)
                bMassCase = false;
            end
            
            %Check that the model has data to write
            %   - 'Children' doesn't count
            %   - Don't write MassGroup data as this has to entered into
            %     the base model bulk data. Difficult to automate!
            tok = obj(1).PropNames;
            tok = tok(~ismember(tok, {'Children', 'MassGroups'}));
            hasData2Write = arrayfun(@(x) checkContents(x, tok), obj);
            
            %Define the part name
            nm = {obj.Name};
            idx = cellfun(@isempty, nm);
            if bCollect
                pName = 'ModelBulkData';
            elseif all(idx)
                pName = '_UnknownPart_';
            else
                %Use the first available name
                ind   = find(~idx, 1, 'first');
                pName = obj(ind).Name;
            end
            
            %Make the file name
            fNam = [pName, '.bdf'];
            
            %Write the data
            if any(hasData2Write)
                
                %Make the file
                fid = fopen(fullfile(dn, fNam), 'w');
                
                %Grab the data, stack it (account for object arrays) &
                %remove any empties
                data = get(obj, tok);
                data = arrayfun(@(i) horzcat(data{:, i}), 1 : size(data, 2), 'Unif', false);
                idx  = cellfun(@isempty, data);
                data(idx) = [];
                tok(idx)  = [];
                cls       = cellfun(@(x) class(x), data, 'Unif', false)';
                
                %Are we writing the design model to a seperate file?
                if bDesModel
                    
                    %Make the file
                    desModelFile = [pName, '_DesignModel.bdf'];
                    dmFid = fopen(fullfile(dn, desModelFile), 'w');
                    obj.writeFileStamp(dmFid);
                    
                    %Grab any designed properties
                    idx = contains(cls, {'awi.fe.opt.PropRelation', ...
                        'awi.fe.opt.PropRelationEqn'});
                    PropRel = horzcat(data{idx});
                    DesModelProps = unique([PropRel.Property]);
                    
                    %Strip these properties from the model
                    prpCls   = arrayfun(@class, DesModelProps, 'Unif', false);
                    ucls     = unique(prpCls);
                    ind      = find(ismember(cls, ucls));
                    desModel = cell(1, numel(ind));
                    for ii = 1 : numel(ind)
                        idx_ = ~ismember(data{ind(ii)},  DesModelProps(ismember(prpCls, ucls{ii})));
                        desModel{ii}  = data{ind(ii)}(~idx_);
                        data{ind(ii)} = data{ind(ii)}(idx_);
                    end
                    desBulkNames = tok(ind);
                    data(idx) = [];
                    tok(idx)  = [];
                    
                    %Grab design variables
                    if bWriteOpt
                        idx = contains(cls, 'awi.fe.DesignVariable');
                        desModel{end + 1}     = data{idx};
                        desBulkNames{end + 1} = tok(idx);
                        data(idx) = [];
                        tok(idx)  = [];
                    end
                    
                    %Write to the file
                    if bComment
                        writePartSummary(dmFid, desModel, desBulkNames, desModelFile);
                    end
                    for ii = 1 : numel(desModel)
                        writeDataToFile(dmFid, desModel{ii}, bComment);
                    end
                    
                    %Close the file
                    fclose(dmFid);
                    
                else
                    
                    desModelFile = [];
                    
                end
                
                %Are we writing the optimisation data?
                if ~bWriteOpt
                    %Find it and remove it
                    idx  = contains(cls, 'awi.fe.opt');
                    data = data(~idx);
                    tok  = tok(~idx);
                end
                
                %Stamp each file with key information
                obj.writeFileStamp(fid);
                
                
                
                %Write a summary of the contents of the part
                if bComment                
                    if ~isempty(obj.Comment)
                        obj.writeComment(obj.Comment,fid)
                    end
                    writePartSummary(fid, data, tok, pName);
                end
                
                %Are we writing mass data?
                if bMassCase
                    fprintf(fid, 'BEGIN massid = %i  label = ''%s''\r\n', ...
                        obj.ID, obj.Name);
                end
                
                %Each property has its own format
                for ii = 1 : numel(data)
                    writeDataToFile(fid, data{ii}, bComment);
                end
                
                %Close the file
                fclose(fid);
                
            else
                
                %Use an empty string to show that the file wasn't generated
                fNam = '';
                desModelFile = [];
                
            end
            
            %Write the child models?
            if ~bCollect
                
                ch = [obj.Children];
                
                %Want to have each child in its own file
                names = cell(size(ch));
                for iCh = 1 : numel(ch)
                    names{iCh} = writeToFile(ch(iCh), dn, bComment, bCollect, bWriteOpt);
                end
                
                %Collect all filenames for use in the header file
                fNam = [{fNam, desModelFile}, horzcat(names{:})];
                
            else
                
                fNam = {fNam, desModelFile};
                
            end
            
            %Strip any empty file names
            fNam = fNam(~cellfun(@isempty, fNam));
            
            function writePartSummary(fid, data, tok, partName)
                %writePartSummary Prints a summary of the contents of the
                %part bulk data into the file with identifier 'fid'.
                %
                % The summary details the number of bulk data entries for
                % each type.
                
                nData   = cellfun(@numel, data, 'Unif', false);
                nChar   = cellfun(@numel, tok);
                maxChar = max(nChar) + 1;
                summary = [tok ; nData];
                
                fprintf(fid, '$\n$  Summary for component %s:\n$\n', partName);
                fprintf(fid, ['$\t* %-', num2str(maxChar), 's : %-8i\n'], summary{:});
                fprintf(fid, '$\n');
                
            end
            
            function writeDataToFile(fid, BulkDataObj, bComment)
                
                
                %Account for heterogeneous arrays
                if strcmp(class(BulkDataObj), 'awi.fe.FEBaseClass') %#ok<STISA>
                    cls_ = arrayfun(@(x) class(x), BulkDataObj, 'Unif', false);
                    uCls = unique(cls_);
                    cellfun(@(x) writeToFile( ...
                        BulkDataObj(ismember(cls_, x)), fid, bComment), uCls);
                else
                    writeToFile(BulkDataObj, fid, bComment);
                end
                
            end
            
        end
    end
    
    methods % visualisation
        function hg = draw(obj, ha, varargin)
            %draw Handles the visualisation of the finite element model.
            %
            % Parameter inputs:
            %   * 'PartList' : List of FE components that are to be drawn.
            %   'PartList' must be a valid 'cellstr' and the names of the
            %   FE components should match the property names. e.g. To draw
            %   all coordinate systems and nodes the user should specify
            %   'PartList', {'CoordSys', 'Nodes'}.
            %
            %   * 'StructureOnly' : If the value of this parameter is
            %   set to 'true' then only the structral part of the FE model
            %   will be drawn.
            %   Default = false
            %
            %   * 'AeroOnly' : If the value of this parameter is set to
            %   'true' then only the aerodynamic parts of the FE model will
            %   be drawn.
            %   Default = false
            %
            % Notes:
            %   * If both 'StructureOnly' and 'AeroOnly' are set to false
            %   then all parts of the FE model will be drawn.
            %   * If both 'StructureOnly' and 'AeroOnly' are set to true
            %   then the structure will override the aero.
            
            hg = [];
            
            %Parse
            p = inputParser;
            addParameter(p, 'PartList'     , []   , @iscellstr);
            addParameter(p, 'StructureOnly', false, @(x)validateattributes(x, {'logical'}, {'scalar'}));
            addParameter(p, 'AeroOnly'     , false, @(x)validateattributes(x, {'logical'}, {'scalar'}));
            parse(p, varargin{:});
            
            %Determine what to draw
            bStructure = true;
            bAero      = true;
            if p.Results.StructureOnly %Draw structure 
                bAero = false;
            end
            if p.Results.AeroOnly      %Draw aero      
                bStructure = false;
            end
            if p.Results.StructureOnly && p.Results.AeroOnly %Ambiguous 
                 bStructure = true;
                 bAero      = false;
            end
            
            %Has the user provided a list of parts?
            if isempty(p.Results.PartList)
                pNames = obj(1).map(:, 2);
            else
                pNames = p.Results.PartList;
                pNames = pNames(ismember(pNames, obj(1).map(:, 2)));
                if isempty(pNames)
                    return
                end
            end
            
            %Create the figure
            if nargin < 2 || isempty(ha)
                hF = figure('Name', 'FE Model');
                ha = axes('Parent', hF, 'NextPlot', 'add');
                axis(ha, 'equal');
            end
            
            %Find the axes
            if isa(ha, 'matlab.ui.Figure')
                ha_ = findobj(ha, 'Type', 'Figure');
                if isempty(ha_)
                    ha = axes('Parent', ha, 'NextPlot', 'add');
                else
                    ha = ha_;
                end
            end
            
            %Validate
            assert(isa(ha, 'matlab.graphics.axis.Axes'), ['Expected ', ...
                'the second argument to ''draw'' to be a valid '     , ...
                '''matlab.graphics.axis.Axes'' object.']);
            
            %Collect all FE models in the hierachy to vectorise
            fem = unique(flatlist(obj));
            
            %'Connections' & the aero panels are treated differently
            idxAero = ismember(pNames, {'Connections', 'AeroPanels', ...
                'AeroControlSurf'});
            if ~any(idxAero)
                bAero = false;
            end
            pNames(idxAero) = [];
            
            %Draw each set of FE objects
            hg = cell(1, numel(pNames));
            if bStructure
                for iP = 1 : numel(pNames)
                    %Gather the data
                    fe = [fem.(pNames{iP})];
                    %Check for empty fe objects
                    if isempty(fe)
                        continue
                    end
                    %Draw it
                    hg{iP} = drawElement(fe, ha);
                end
                %Joints & Bush Elements
                if ~isempty([fem.Connections])
                    con  = [fem.Connections];
                    idxJ = arrayfun(@(x) isa(x, 'awi.fe.Joint'), con);
                    idxB = arrayfun(@(x) isa(x, 'awi.fe.BushElement'), con);
                    hJ   = drawElement(con(idxJ), ha); %Draw the joints
                    hBu  = drawElement(con(idxB), ha); %Draw the bush elements
                    hg{end + 1} = [hJ ; hBu];
                end
            end
            
            %Aerodynamic Panels & Control Surfaces
            %   - Grab all the panel objects and check if any of them have
            %     been defined as control surfaces
            panels        = [fem.AeroPanels];
            control       = [fem.AeroControlSurf];
            if bAero && ~and(isempty(control), isempty(panels)) 
                if isempty(control) && ~isempty(panels)
                    hg{end + 1} = drawElement(panels, ha);
                else
                    panelSet      = [control.AeroPanelSet];
                    controlPanels = [panelSet.AeroPanels];
                    idxControl    = ismember(panels, controlPanels);
                    hg{end + 1}   = drawElement(panels(~idxControl), ha);
                    hg{end + 1}   = drawElement(controlPanels, ha  , ...
                        'Tag'      , 'Aerodynamic Control Surfaces', ...
                        'FaceColor', ([249, 189, 24] ./ 255));
                end
            end
            
            %Collapse cell-array
            hg = vertcat(hg{:});
            
            %Do not allow FEM to be highlighted
            set(hg, 'SelectionHighlight', 'off');
            
            %Update legend
            legend(ha, hg, get(hg, 'Tag'), ...
                'ItemHitFcn', @i_toggle_visible);
            
            function i_toggle_visible(~, evt)
                
                switch evt.Peer.Visible
                    
                    case 'on'
                        evt.Peer.Visible = 'off';
                    case 'off'
                        evt.Peer.Visible = 'on';
                end
                
            end
            
        end
        function hg = drawMassDistribution(obj, ha, varargin)
            %drawMassDistribution Draws the mass distribution of the FEM.
            %
            % Each point mass object is represented as a sphere whose
            % radius is scaled to with respect to the largest mass
            
            hg = [];
            
            %Create the figure
            if nargin < 2 || isempty(ha)
                hF = figure('Name', 'FE Model');
                ha = axes('Parent', hF, 'NextPlot', 'add');
                axis(ha, 'equal');
            end
            
            %Find the axes
            if isa(ha, 'matlab.ui.Figure')
                ha_ = findobj(ha, 'Type', 'Figure');
                if isempty(ha_)
                    ha = axes('Parent', ha, 'NextPlot', 'add');
                else
                    ha = ha_;
                end
            end
            
            %Validate
            assert(isa(ha, 'matlab.graphics.axis.Axes'), ['Expected ', ...
                'the second argument to ''draw'' to be a valid '     , ...
                '''matlab.graphics.axis.Axes'' object.']);
            
            %Default colour
            clr = [21, 176 237]./255;
            
            %Parse
            p = inputParser;
            addParameter(p, 'MassColor'     , clr , @(x)validateattributes(x, {'numeric'}, {'nonnegative', 'row', 'ncols', 3, '<=', 1}));
            addParameter(p, 'ScaleMasses'   , true, @(x)validateattributes(x, {'logical'}, {'scalar'}));
            addParameter(p, 'MaxRadius'     , 1   , @(x)validateattributes(x, {'numeric'}, {'scalar', 'positive'}));
            addParameter(p, 'MassCases'     , []  , @(x)validateattributes(x, {'awi.fe.MassModel'}, {'vector'}));
            addParameter(p, 'MassCaseNames' , []  , @iscellstr);
            addParameter(p, 'MassCaseColors', []  , @(x)validateattributes(x, {'numeric'}, {'nonnegative', '2d', 'ncols', 3, '<=', 1}));
            parse(p, varargin{:});
            assert(numel(p.Results.MassCaseNames) == numel(p.Results.MassCases), ...
                ['The number of user-supplied "MassCaseNames" should ', ...
                'equal the number of user-supplied "MassCases".']);
            assert(numel(p.Results.MassCaseNames) == numel(p.Results.MassCases), ...
                ['The number of user-supplied "MassCaseColors" should ', ...
                'equal the number of user-supplied "MassCases".']);
            
            %Collect all FE models in the hierachy to vectorise
            fem = unique(flatlist(obj));
            
            %Check for point masses
            PointMasses = [fem.PointMasses];
            if isempty(PointMasses)
                warning(['No ''awi.fe.PointMass'' objects in the ', ...
                    'finite element model. Unable to show mass '  , ...
                    'distribution']);
                return
            end
            
            %Continue with point masses that have mass and position
            idxM    = arrayfun(@(m) ~isempty(m.M)    , PointMasses);
            idxNode = arrayfun(@(m) ~isempty(m.Node) , PointMasses);
            PointMasses = PointMasses(and(idxM, idxNode));
            idxPos  = arrayfun(@(m) ~isempty(m.Node.X), PointMasses);
            PointMasses = PointMasses(idxPos);
            if isempty(PointMasses)
                return
            end
            
            %Any beam objects to show skeleton of structure?
            Beams = [fem.Beams];
            if ~isempty(Beams)
                hg = drawElement(Beams, ha);
            end
            
            %Draw the masses as spheres
            hg(end + 1) = i_drawMassAsSpheres(ha, PointMasses, p.Results.ScaleMasses, p.Results.MaxRadius);
            set(hg(end), 'FaceColor', p.Results.MassColor);
            
            %Any MassCases?
            MassCases = (p.Results.MassCases);
            if ~isempty(MassCases)
                nMC = numel(MassCases);
                nam = p.Results.MassCaseNames;
                if isempty(nam)
                    nam = arrayfun(@(i) sprintf('Mass Case %i', 1 : nMassCase, 'Unif', false));
                end
                clrs = p.Results.MassCaseColors;
                hg_ = gobjects(1, nMC);
                %Recurse throguh validation and then plotting...
                for i = 1 : nMC
                    hg_(i) = drawMassDistribution(MassCases(i), ha, ...
                        'MassColor'  , clrs(i, :), ...
                        'MaxRadius'  ,  p.Results.MaxRadius, ...
                        'ScaleMasses',  p.Results.ScaleMasses);
                end
                hg = [hg, hg_];
            end
            
            %Do not allow FEM to be highlighted
            set(hg, 'SelectionHighlight', 'off');
            
            function hg = i_drawMassAsSpheres(ha, PointMasses, bScaleMasses, maxR)
                %i_drawMassAsSpheres Draws each 'awi.fe.PointMass' object
                %as a sphere with a radius that is (optionally) scaled with
                %respect to its mass.
                
                %Grab mass data
                r     = get([PointMasses.Node], {'X'});
                rOff  = [PointMasses.X];
                r     = horzcat(r{:}) + rOff;
                m     = [PointMasses.M];
                nMass = numel(m);
                
                %Determine readius of spheres
                if bScaleMasses
                    rad = maxR .* m ./ max(m);
                else
                    rad = repmat(maxR, [1, nMass]);
                end
                
                %Coordinates of a unit sphere at (0, 0, 0)
                [x,y,z] = sphere;
                fvc     = surf2patch(x, y, z);
                
                %Generate single set of patch coordinates for vectorised plot
                
                %   - Vertex index for each sphere
                nVtx  = size(fvc.vertices, 1);
                iVtx  = [0, repmat(nVtx, [1, (nMass - 1)])];
                iVtx  = cumsum(iVtx);
                %   - Repeat the face and vertex data
                vtx = repmat(fvc.vertices, [1, 1, nMass]);
                fcs = repmat(fvc.faces   , [1, 1, nMass]);
                %   - Scale radius and add origin
                vtx = vtx .* permute(rad, [1, 3, 2]);
                vtx = vtx + permute(r, [3, 1, 2]);
                %   - Set correct vertex index number
                fcs = fcs + permute(iVtx, [3, 1, 2]);
                %   - Stack coordinates ready for plotting
                vtx = arrayfun(@(i) vtx(:, :, i), 1 : size(vtx, 3), 'Unif', false);
                fcs = arrayfun(@(i) fcs(:, :, i), 1 : size(fcs, 3), 'Unif', false);
                vtx = vertcat(vtx{:});
                fcs = vertcat(fcs{:});
                
                %Plot
                hg = patch(ha, 'Faces', fcs, 'Vertices', vtx);
                
            end
            
        end
        function hg = drawBeamOffsets(obj, ha)
            hg = [];
            
            %Create the figure
            if nargin < 2 || isempty(ha)
                hF = figure('Name', 'FE Model');
                ha = axes('Parent', hF, 'NextPlot', 'add');
                axis(ha, 'equal');
            end
            
            %Find the axes
            if isa(ha, 'matlab.ui.Figure')
                ha_ = findobj(ha, 'Type', 'Figure');
                if isempty(ha_)
                    ha = axes('Parent', ha, 'NextPlot', 'add');
                else
                    ha = ha_;
                end
            end
            
            %Validate
            assert(isa(ha, 'matlab.graphics.axis.Axes'), ['Expected ', ...
                'the second argument to ''draw'' to be a valid '     , ...
                '''matlab.graphics.axis.Axes'' object.']);
            
            %Parse
            %             p = inputParser;
            %             addParameter(p, 'ScaleMasses'   , true, @(x)validateattributes(x, {'logical'}, {'scalar'}));
            %             addParameter(p, 'MaxRadius'     , 1   , @(x)validateattributes(x, {'numeric'}, {'scalar', 'positive'}));
            %             addParameter(p, 'MassCases'     , []  , @(x)validateattributes(x, {'awi.fe.MassModel'}, {'vector'}));
            %             addParameter(p, 'MassCaseNames' , []  , @iscellstr);
            %             addParameter(p, 'MassCaseColors', []  , @(x)validateattributes(x, {'numeric'}, {'nonnegative', '2d', 'ncols', 3, '<', 1}));
            %             parse(p, varargin{:});
            %             assert(numel(p.Results.MassCaseNames) == numel(p.Results.MassCases), ...
            %                 ['The number of user-supplied "MassCaseNames" should ', ...
            %                 'equal the number of user-supplied "MassCases".']);
            %             assert(numel(p.Results.MassCaseNames) == numel(p.Results.MassCases), ...
            %                 ['The number of user-supplied "MassCaseColors" should ', ...
            %                 'equal the number of user-supplied "MassCases".']);
            %
            %Collect all FE models in the hierachy to vectorise
            fem = unique(flatlist(obj));
            
            %Check for point masses
            Beams = [fem.Beams];
            if isempty(Beams)
                warning(['No ''awi.fe.Beam'' objects in the finite ', ...
                    'element model. Unable to show offset data.']);
                return
            end
            
            %Continue with point masses that have mass and position
            idxNode = arrayfun(@(b) ~isempty(b.Nodes), Beams);
            idxX    = arrayfun(@(b) ~isempty(b.X)    , Beams);
            Beams   = Beams(and(idxX, idxNode));
            idxPos  = arrayfun(@(b) any(arrayfun(@(n) ~isempty(n.X), b.Nodes)), Beams);
            Beams = Beams(idxPos);
            if isempty(Beams)
                return
            end
            
            %Grab coordinates of end-A
            Nodes = [Beams.Nodes];
            rA = [Nodes(1, :).X];
            
            %Calculate orientation vectors of beam local coordinate system
            dR = arrayfun(@(b) diff([b.Nodes.X], [], 2), Beams, 'Unif', false);
            dR = horzcat(dR{:});
            v1 = dR ./ sqrt((dR(1, :).^2 + dR(2, :).^2 + dR(3, :).^2));
            v2 = [Beams.X];
            
            %Draw beam objects
            hg = drawElement(Beams, ha);
            
            %Draw orientation vectors
            xyz1 = rA + v1;
            x1 = [rA(1, :) ; xyz1(1, :)];
            y1 = [rA(2, :) ; xyz1(2, :)];
            z1 = [rA(3, :) ; xyz1(3, :)];
            plot3(ha, x1, y1, z1, 'r-', 'LineWidth', 2);
            %
            xyz2 = rA + v2;
            x2 = [rA(1, :) ; xyz2(1, :)];
            y2 = [rA(2, :) ; xyz2(2, :)];
            z2 = [rA(3, :) ; xyz2(3, :)];
            plot3(ha, x2, y2, z2, 'g-', 'LineWidth', 2);
            
            %Draw the masses as spheres
            hg(end + 1) = i_drawMassAsSpheres(ha, Beams, p.Results.ScaleMasses, p.Results.MaxRadius);
            set(hg(end), 'FaceColor', [21, 176 237]./255);
            
            %Any MassCases?
            MassCases = (p.Results.MassCases);
            if ~isempty(MassCases)
                nMC = numel(MassCases);
                nam = p.Results.MassCaseNames;
                if isempty(nam)
                    nam = arrayfun(@(i) sprintf('Mass Case %i', 1 : nMassCase, 'Unif', false));
                end
                %Recurse throguh validation and then plotting...
                %                 hg_ = arrayfun(@(mc)
            end
            
            %Do not allow FEM to be highlighted
            set(hg, 'SelectionHighlight', 'off');
            
        end
        function hg = drawResults(obj, ha, ResultsObj, sf)
            %drawResults Draws the results data contained in 'ResultsObj'
            %in the axis 'ha' for the FEM 'obj'.
            
            hg = [];
            
            if nargin < 4
                sf = 1;
            end
            
            %Create the figure
            if nargin < 2 || isempty(ha)
                hF = figure('Name', 'FE Model');
                ha = axes('Parent', hF, 'NextPlot', 'add');
                axis(ha, 'equal');
            end
            
            %Find the axes
            if isa(ha, 'matlab.ui.Figure')
                ha_ = findobj(ha, 'Type', 'Figure');
                if isempty(ha_)
                    ha = axes('Parent', ha, 'NextPlot', 'add');
                else
                    ha = ha_;
                end
            end
            
            %Validate
            assert(isa(ha, 'matlab.graphics.axis.Axes'), ['Expected ', ...
                'the second argument to ''draw'' to be a valid '     , ...
                '''matlab.graphics.axis.Axes'' object.']);
                        
            %Attach results 
            %Assign to the Node objects
            set(ResultsObj.Nodes, {'GlobalTranslation'}, num2cell(ResultsObj.Translation(:, :, end) * sf, 1)');
            set(ResultsObj.Nodes, 'DrawMode', 'deformed');
            hg = draw(obj, ha);
            
            %Detach results
            
        end
    end
    
    methods % analytical methods
        function [beamEta, beamR] = calculateBeamEta(obj)
            %calcualteBeamEta Calculates the piecewise-linear 'eta' value
            %(i.e. along the beam) of various beam elements in the model
            %asssuming that they are connected end-to-end in order.
            
            assert(numel(obj) == 1, ['Method ''calculateBeamEta'' ', ...
                'is not valid for object arrays.']);
            
            beamEta = [];
            
            %Grab nodes
            Nodes = obj.BeamNodes; %#ok<*PROP>
            
            if isempty(Nodes)
                return
            end
            
            %Assume beams are arranged in order along the beam
            Nodes = [Nodes(1, :), Nodes(2, end)];
            
            if any([Nodes.CP] ~= 0)
                error('Update for nodes defined in the non-basic coordinate system.');
            end
            
            %Calculate length of line through all beam nodes
            X = [Nodes.X];
            beamR = awi.model.Stick.getLineLength(X(1, :), X(2, :), X(3, :));
            
            %Normalise
            beamEta = beamR ./ beamR(end);
            
            %End-A & End-B
            beamEta = [beamEta(1 : end - 1) ; beamEta(2 : end)];
            
        end
        function assignMassNodes(obj, PointMass, massPos)
            %assignMassNodes Assigns the 'Node' property of the 'PointMass'
            %objects based on nearest 'BeamNode' to that point mass.
            
            validateattributes(PointMass, {'awi.fe.PointMass'}, {'row'}, ...
                class(obj), 'PointMass');
            validateattributes(massPos, {'numeric'}, {'2d', 'ncols', 3}, ...
                class(obj), 'massPos');
            
            %Position of beam nodes
            BeamNodes = obj.BeamNodes; %#ok<*PROPLC>
            if isempty(BeamNodes) %Escape route
                return
            end
            BeamNodes = unique(BeamNodes(:), 'stable')';
            %             beamX = [BeamNodes.X];
            %             x = beamX(1, :);
            %             y = beamX(2, :);
            %             z = beamX(3, :);
            
            %Position of masses
            xm = massPos(:, 1);
            ym = massPos(:, 2);
            zm = massPos(:, 3);
            
            %Find closest nodes
            NearNode = obj.findClosestNode(BeamNodes, massPos');
            
            %             %Find the nearest node to the mass coordinates and assign
            %             %that as the 'Node' property.
            %             r    = sqrt((x - xm).^2 + (y - ym).^2 + (z - zm).^2);
            %             rMin = min(r, [], 2);
            %             index    = arrayfun(@(i) find(r(i, :) == rMin(i), 1), 1 : numel(rMin));
            %             NearNode = BeamNodes(index);
            
            
            %Assign Nodes for each mass
            set(PointMass, {'Node'}, num2cell(NearNode)');
            
            %Work out the offset distance
            nodeCoords = [NearNode.X]';
            dx = xm - nodeCoords(:, 1);
            dy = ym - nodeCoords(:, 2);
            dz = zm - nodeCoords(:, 3);
            Xoff = [dx, dy, dz]';
            set(PointMass, {'X'}, num2cell(Xoff, 1)');
            
        end
    end
    
    methods % writing progress
        function log(obj, str)
            %log Sends the message 'str' to the log function and prefaces
            %it with the current indentation character.
            msg = [obj.IndentChar, str];
            obj.LogFcn(msg);
        end
    end
    
    methods (Static) %Helper functions
        function NearNode = findClosestNode(Nodes, coords)
            %findClosestNode Finds the closest 'awi.fe.Node' objects to the
            %coordinates in 'coords'.
            
            NearNode = [];
            
            %Parse
            validateattributes(Nodes , {'awi.fe.Node'}, {'row'}           , 'findClosestNode', 'Nodes');
            validateattributes(coords, {'numeric'}    , {'2d', 'nrows', 3}, 'findClosestNode', 'coords');
            
            %Check for empty coordinate data
            nodeX = get(Nodes, {'X'});
            idx   = ~cellfun(@isempty, nodeX);
            Nodes = Nodes(idx);
            if isempty(Nodes)
                return
            end
            
            %Get coordinates of nodes
            nodeX = horzcat(nodeX{:});
            xn = nodeX(1, :);
            yn = nodeX(2, :);
            zn = nodeX(3, :);
            
            %Transpose coordinate matrix so we can compare position of all
            %node coordinates with the requested coordinates.
            coords = coords';
            
            %Position of masses
            xc = coords(:, 1);
            yc = coords(:, 2);
            zc = coords(:, 3);
            
            %Find the nearest node to the mass coordinates and assign
            %that as the 'Node' property.
            r    = sqrt((xn - xc).^2 + (yn - yc).^2 + (zn - zc).^2);
            rMin = min(r, [], 2);
            index    = arrayfun(@(i) find(r(i, :) == rMin(i), 1), 1 : numel(rMin));
            NearNode = Nodes(index);
            
        end
        function [Nodes, Joint, Bush, BushPrp] = makeDefaultConnection
            %makeDefaultConnection Creates the default FE objects to
            %faciliate a semi-rigid link between two different FE models.
            %
            % Default is :
            %   - 2 x Nodes
            %   - 2 x Joint
            %   - 1 x BushElement
            %   - 1 x BushProp
            
            Nodes   = arrayfun(@(~) awi.fe.Node , 1 : 2);
            Joint   = arrayfun(@(~) awi.fe.Joint, 1 : 2);
            Bush    = awi.fe.BushElement;
            BushPrp = awi.fe.BushProp;
            
            %Set bush stiffness to be quite high
            k    = 1e11;
            kNam = arrayfun(@(i) ['K', num2str(i)], 1 : 6, 'Unif', false);
            kVal = repmat({k}, [1, 6]);
            set(BushPrp, kNam, kVal);
            
        end
    end
    
    methods (Access = protected) % custom display
        function s = getHeader(obj)
            %getHeader Returns the title header for the object display.
            
            s   = getHeader@matlab.mixin.CustomDisplay(obj);
            ind = strfind(s, 'with properties');
            if ~isempty(ind)
                s   = [s(1 : ind - 1), 'with non-empty FE properties:'];
            end
            
        end
        function groups = getPropertyGroups(obj)
            %getPropertyGroups Returns a 'matlab.mixin.util.PropertyGroup'
            %object for displaying the object's properties.
            %
            % Detailed Description:
            %   - Filters the list of FE properties to return only those
            %     that are populated.
            
            %Initiate the group
            groups = matlab.mixin.util.PropertyGroup;
            
            %What properties have we got?
            allPrp = properties(obj);
            fePrp  = obj(1).map(:, 2);
            
            %Always display non-FE properties
            prp    = allPrp(~ismember(allPrp, fePrp));            
            
            %What is the state of the FE properties?
            fePrpVal = get(obj, fePrp);
            
            %Only keep non-empty properties
            bNonEmpty = any(~cellfun(@isempty, fePrpVal), 1);            
            fePrp     = fePrp(bNonEmpty);
            fePrpVal  = fePrpVal(:, bNonEmpty);
            
            %Consolidate
            prps = [prp ; fePrp];
            [prps, ia, ~] = unique(prps, 'stable');
            if numel(obj) > 1
                groups.PropertyList = prps;
            else                
                prpVal = [get(obj, prp) , fePrpVal];                
                prpVal         = prpVal(ia);
                groups.PropertyList = cell2struct(prpVal, prps', 2);
            end
                        
        end
    end
    
end

