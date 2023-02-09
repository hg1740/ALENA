classdef (ConstructOnLoad) Beamable < matlab.mixin.SetGet & mvc.mixin.Dynamicable 
    %Beamable Defines a generic beam which can have properties along its
    %length.
    %
    % The 'Beamable' class is dependent on the 'BeamProperty' and
    % 'BeamPropertyObject' classes and will not function correctly without
    % it.
    
    %Helper properties for storing 'awi.mixin.BeamProperty' objects as well
    %as generic objects that are located along the beam.
    properties (Hidden, SetAccess = private)
        %Array of 'awi.model.BeamProperty' objects which are defined along
        %the span of the beam
        BeamProperties
        %Name of the beam objects that have been assigned
        BeamObjectNames = {}; 
        %Classes of the beam objects that have been assigned - for
        %validation
        BeamObjectClass = {};
        %Name of the properties that are dependent on the beam objects
        BeamObjectProps = {};
    end
    
    %List of beam property types
    properties (Dependent)
        %List of all beam properties assigned to this object
        BeamPropertyNames
        %List of all unique beam property types
        BeamPropertyTypes        
    end
    
    %Dynamic Property Handles and Listeners - TODO : Eventually, we should
    %replace this with the "listeners('set')" command and allow Collectable
    %to handle everything
    properties (Hidden, SetAccess = private)
       hDynListener
       hDynBeamProp
    end
    
    %Helper properties for defining the names of dynamic properties and
    %their corresponding properties in the 'BeamProperty', class.
    properties (Hidden, Constant, Access = private)
        %Each beam property that is added to the object will be appended by
        %the suffixes defined in 'DynPropTag'. For each beam property that 
        %is added to the model a total of numel(DynPropTag) properties are
        %added to the object as dynamic properties.
        DynPropTag = {'', '_eta', '_var', '_flag', '_'}; %TODO - Change '_' to '_int'
        %'DynPropKey' defines the property name in the 'BeamProperty' class
        %that corresponds to the tag in 'DynPropTag'.
        DynPropKey = {'Value', 'Distribution', 'Variation', 'AxisFlag', 'Interpolated'};
        %'DynPropObs' controls which of the dynamic properties are
        %observable.
        DynPropObs = [true, true, true, true, true];
        %'DynPropHid' controls which of the dynamic properties are hidden
        DynPropHid = [false, false, false, false, true];
    end
    
    %Helper properties for defining the names of dynamic beam object
    %properties
    properties (Hidden, Constant, Access = private)
        %Each beam property object that is added to the object will be
        %appended by the suffixes defined in ''BeamObjTag''. For each beam
        %property object that is added to the model a total of
        %numel(BeamObjTag) properties are added to the object as dynamic
        %properties.
        BeamObjTag = {'', '_eta', '_flag', '_var'}
        %'DynPropObs' controls which of the dynamic beam property objects 
        %are observable.
        BeamObjObs = [true, true, true, true];
    end
      
    methods % set / get
        function set.BeamProperties(obj, val)     %set.BeamProperties    
            %set.BeamProperties Set method for the property 'BeamProperties'
            %
            % 'BeamProperties' must be a row vector of class
            % 'awi.mixin.BeamProperty'.
            validateattributes(val, {'awi.mixin.BeamProperty'}, {'row'}, ...
                class(obj), 'BeamProperties');
            obj.BeamProperties = val;
        end  
        function val = get.BeamPropertyNames(obj) %get.BeamPropertyNames 
            val = [];
            BP = obj.BeamProperties;
            if isempty(BP)
                return
            end
            val = {BP.Name};
        end
        function val = get.BeamPropertyTypes(obj) %get.BeamPropertyTypes 
            %get.BeamPropertyTypes Get method for the dependent property
            %'BeamPropertyTypes'.
            %
            %   - 'BeamPropertyTypes' is a cellstr of all unique beam
            %      property types in the object.
            
            if isempty(obj.BeamProperties) %Escape route
                val = [];
                return
            end
            
            %Simple
            allTypes = {obj.BeamProperties.Type};
            val      = unique(allTypes);
        end
    end
    
    methods % constructor
        function obj = Beamable
            %Beamable Constructor method for the class 'Beamable'.
            
            %Add a link to the beam geometry
            if isa(obj, 'awi.model.Stick')
%                 listeners(obj, 'set', @obj.updateBeamGeometry
%                 addListeners(obj, 'XData', @obj.updateBeamGeometry, 'PostSet');
                %Listen for changes in the beam geometry
                addlistener(obj, {'XData', 'YData', 'ZData'}, ...
                    'PostSet', @obj.updateBeamGeometry);
            end
           
            %Add a link to the 'SpanVector' property
            if isa(obj, 'awi.model.LiftingSurface')
                %Listen for changes in the span vector
                addlistener(obj, 'SpanVector', ...
                    'PostSet', @obj.updateBeamPropertyAxis);
            end
        end
    end
    
    methods (Sealed) % interacting with beam properties 
        % defining beam properties
        function addBeamProperty(obj, name, varargin)
            %addBeamProperty Adds an instance of 'awi.mixin.BeamProperty'
            %to the 'BeamProperties' list in the class 'obj' and
            %initialises that 'Quantity' property with the value of 'name'.
            %
            % In order to facilite a smooth user-interface with the object
            % this method will initialise a group of 'public-facing'
            % dynamic properties that provide a shortcut to the true beam
            % properties which themselves are hidden.
            %
            % For this reason, the implementation of 'Beamable' is tightly
            % linked with the 'awi.mixin.BeamProperty' class and will not
            % work without it.
            %
            %
            % TODO - Add additional arguments for accepting beam properties
            % that are not row vectors. i.e. Shear offsets in (y,z)
            % directons and stiffness matrices that are 6x6 etc. What about
            % objects such as the CoordSys class?
            %
            
            %Create an instance of 'BeamProperty' & initialise the name
            BP = awi.mixin.BeamProperty('Quantity', name, varargin{:});
            
            %Add to the list of beam properties
            obj.BeamProperties = [obj.BeamProperties, BP];
            
            %Create a series dynamic property with the same name that acts 
            %as the visible beam property.
            p = makeDynProps(obj, name, get(BP, obj.DynPropKey));
            
            %We want to listen to changes in this property so we can relate
            %back to the handle array of 'BeamProperty' classes.
            %   - Set & Get Observable
            %   - Abort Set
            %   - Hidden *Only applies to some of the properties.
            %   - Dependent *Only the 'interpolated' property is dependent.
            idx  = ismember(obj.DynPropKey, 'Interpolated');
            temp = num2cell(obj.DynPropObs);
            [p.SetObservable] = deal(temp{:});
            [p.GetObservable] = deal(temp{:});
            [p.AbortSet]      = deal(temp{:});
            temp = num2cell(obj.DynPropHid);
            [p.Hidden]        = deal(temp{:});
            p(idx).Dependent  = true;
            clear temp
            
            %For members of a collection
            if isa(obj, 'mvc.mixin.Collectable')
                %TODO - Make this work!!! (Chris S, code properly you
                %douche!)
                %                 listeners(obj, 'set');
            end
            
            %Define listener for 'PostSet' on each new dynamic property
            %   - Only do this for those objects that are observable!
            lh(1, 1) = addlistener(obj, {p(and(obj.DynPropObs, ~idx)).Name}, 'PostSet', @obj.updateBeamProperty);
            
            %Define listener to retrieve the interpolated value of the beam
            %property whenever the '_' property is requested
            lh(1, 2) = addlistener(obj, [p(1).Name, '_'], 'PreGet', @obj.retrieveBeamProperty);
            
            %             %Make a note
            %             if isa(obj, 'mvc.mixin.Dynamicable')
            %                 %Add to the 'hDynamicProperties' list
            %                 obj.hDynamicProperties = [obj.hDynamicProperties, p];
            %             else
            %Use custom list
            if isempty(obj.hDynBeamProp)
                obj.hDynBeamProp = p;
                obj.hDynListener = lh; %TODO - Decide if I want to keep this. It is not being used currently.
            else
                obj.hDynBeamProp = [obj.hDynBeamProp, p];
                obj.hDynListener = [obj.hDynListener, lh];
            end
            %             end
            %
            %             %Update public facing property using the default values in BP
            %             var = cellfun(@(x) [name, x], obj.DynPropTag, 'Unif', false);
            %             val = get(BP, obj.DynPropKey);
            %             set(obj, var, val);
            
            %Create/Add a property group for the dynamic properties
            obj.addPropertyGroup(BP.Type, ...
                BP.Quantity, BP.Name, BP.Description, {});
            %Add new parameter set to obj.PropertyGroups
            %             tit = obj.ParameterSets(tdx).DisplayName;
            %             c = obj.ParameterSets(tdx).Parameters';
            %             obj.addPropertyGroup(tit, c{:});
            
        end        
        % defining beam property objects
        function addBeamPropObj(obj, cls, var, objName, varargin)
            %addBeamPropObj Adds a set of dynamic properties to the
            %'Beamable' object that allow beam quantities to be defined
            %over the length of the beam based on the data in the beam
            %property of class 'cls'.
            
            %Parse inputs
            if nargin < 3 || isempty(varargin) %Escape route
                %No point adding the dynamic properties if the user hasn't
                %provided enough information!
                return
            end
            assert(ischar(cls), 'Expected the object class name to be a character array.');
            bpo = str2func(cls);
            bpo = bpo();
            try
            assert(isa(bpo, 'awi.model.Entity'), 'Expected the beam property object to be of type ''awi.model.Entity''.');
            catch ME
                disp('Wait')
            end
            awi.mixin.BeamProperty.validateVariation(obj, var);
            if nargin < 4 %Just use the object type if no name is supplied
                objName = obj.Type;
            end           
            assert(iscellstr(varargin), 'Expected the extra property names to be a cell-array of strings.'); 
            
            %User may specify beam property 'Type' as part of input
            idx = ismember(varargin, 'Type');
            if any(idx)
                ind = find(idx);
                type = varargin{ind + 1};
                varargin([ind, ind + 1]) = [];
            else
                type = objName;
            end
            
            %Does the object already have a dynamic property of this name?
            %   TODO : Finish this off
            if ~isempty(obj.BeamObjectClass)
                idx = ismember(obj.BeamObjectClass, cls);
            end
            
            %Create a dynamic property with the same name that acts as the
            %visible beam property
            val = cell(size(obj.BeamObjTag));
            val{ismember(obj.BeamObjTag, '_var')} = var;
            p = arrayfun(@(i)addDynamic(obj, [objName, obj.BeamObjTag{i}], ...
                val{i}), 1 : numel(obj.BeamObjTag), 'Unif', false);
            p = horzcat(p{:});
            
            %We want to listen to changes in these properties so we can 
            %validate the values that the user sets
            %   - Set & Get Observable
            %   - Abort Set
            %   - SetAccess = 'private'
            temp = num2cell(obj.BeamObjObs);
%             sa   = repmat({'protected'}, size(obj.BeamObjObs));
%             sa   = repmat({'private'}, size(obj.BeamObjObs));
            [p.SetObservable] = deal(temp{:});
            [p.GetObservable] = deal(temp{:});
            [p.AbortSet]      = deal(temp{:});
%             [p.SetAccess]     = deal(sa{:});
            
            %Update the list of beam object names and class
            obj.BeamObjectNames{end + 1} = objName;
            obj.BeamObjectClass{end + 1} = cls;
            obj.BeamObjectProps{end + 1} = varargin;
            
            %Add 'PostSet' methods to allow the dynamic properties to be
            %validated
            lh(1, 1) = addlistener(obj, p(1).Name, 'PostSet', @obj.validateBeamPropObj);
            lh(1, 2) = addlistener(obj, {p(2:end).Name}, 'PostSet', @obj.updateSpawnedBeamProperty);
            
            %Use custom list
            if isempty(obj.hDynBeamProp)
                obj.hDynBeamProp = p;
                obj.hDynListener = lh; %TODO - Decide if I want to keep this. It is not being used currently.
            else
                obj.hDynBeamProp = [obj.hDynBeamProp, p];
                obj.hDynListener = [obj.hDynListener, lh];
            end
            
            
            %Make the 'awi.mixin.BeamProperty' objects & set the name
            BP = arrayfun(@(i) awi.mixin.BeamProperty( ...
                'Quantity', varargin{i}, 'Type', type), ...
                1 : numel(varargin), 'Unif', false);
            BP = horzcat(BP{:});
            
            %Override the variation to match the beam property object
            set(BP, 'Variation', var);
            
            %Add to the list of beam properties
            obj.BeamProperties = [obj.BeamProperties, BP];
            
            %Create the dynamic properties which will depend on the beam
            %property object
            %   - Initialise using defaults from the beam properties
            %   - 
            val = get(BP, obj.DynPropKey);
            p = arrayfun(@(i) makeDynProps(obj, varargin{i}, val(i, :)), ...
                1 : numel(varargin), 'Unif', false);
            p = horzcat(p{:});
            
            %We don't want the user to be able to define the value of these
            %properties
            %   - GetObservable
            %   - SetAccess = 'protected'
            temp = repmat({true}, size(p));
%             sa   = repmat({'protected'}, size(p));
            [p.GetObservable] = deal(temp{:});
%             [p.SetAccess]     = deal(sa{:});       
            
            %Define listeners for the 'spawned' beam properties
            nProp = numel(obj.DynPropKey); %No. of spawned beam properties
            %   - 'Value' properties need a 'PreGet' listener which
            %      interrogates the relevant beam property objects and
            %      retrieves the value stored in those objects.
            pValue = p(1 : nProp : end);
%             get_fn = repmat({@obj.deriveBeamProp}, size(pValue));
%             [pValue.GetMethod] = deal(get_fn{:});
            lh(1, 1) = addlistener(obj, {pValue.Name}, 'PreGet', @obj.deriveBeamProperty); 
            %   - 'Distribution' & 'AxisFlag' properties need a 'PostSet' 
            %     listener to pass the value from the beam property object. 
            idx = ismember(obj.DynPropKey, {'Distribution', 'Variation', 'AxisFlag'});
            ind = find(idx);
            ind = arrayfun(@(i) (ind(i) : nProp : numel(p)), 1 : numel(ind), 'Unif', false);
            ind = horzcat(ind{:});
            pDAF = p(ind);
            temp = repmat({true}, size(pDAF));
            [pDAF.SetObservable] = deal(temp{:});
            lh(1, 2) = addlistener(obj, {pDAF.Name}, 'PostSet', @obj.updateBeamProperty);
            
            %Define listener to retrieve the interpolated value of the beam
            %property whenever the '_' property is requested
            pNam = cellfun(@(x) [x, '_'], {pValue.Name}, 'Unif', false);
            lh(1, 3) = addlistener(obj, pNam, 'PreGet', @obj.retrieveBeamProperty);
            
            
             %Use custom list
            if isempty(obj.hDynBeamProp)
                obj.hDynBeamProp = p;
%                 obj.hDynListener = lh; %TODO - Decide if I want to keep this. It is not being used currently.
            else
                obj.hDynBeamProp = [obj.hDynBeamProp, p];
%                 obj.hDynListener = [obj.hDynListener, lh];
            end
            
            
        end        
        % adding/removing generic objects along the span of the beam
        function assignBeamObject(obj, beamObj, eta, mode)
            %assignBeamObject Assigns a beam object to a particular
            %non-dimensional position along the length of the beam.
            %
            % If 'mode' is set to 'replace' then the current distribution
            % is replaced by the new distribution. N.B. The default
            % behaviour is to augment the current distribution.
            
            if nargin < 3 %Escape route 
                return
            end
            
            if nargin < 4 %Default behaviour is to augment 
                mode = 'augment';
            end
            
            %Parse the inputs
            [nam, eNam] = parseBeamPropInputs(obj, beamObj, eta);
            
            if isempty(nam) || isempty(eNam) %Escape route
                return
            end
            
            switch mode
                case 'augment'
                    %Grab current distribution of beam property objects
                    %beamObj = arrayfun(@(o) [o.(nam), beamObj], obj, 'Unif', false);
                    beamObj = [obj.(nam), beamObj];
                    objEta  = [obj.(eNam), eta];
                case 'replace'
                    %Remove the current distribution
                    delete(obj.(nam));
                    obj.(nam) = [];
                    %Just use the new distribution and new objects
                    objEta    = eta;
                otherwise
                    error(['Unknown option %s for the method ', ...
                        '''assignBeamObject'' of object %s']  , ...
                        mode, class(obj));
            end             
                       
            %Sort by eta position
            [objEta, index] = sort(objEta, 'ascend');
            
            %Assign the beam object 
            obj.(nam)  = beamObj(index);
            obj.(eNam) = objEta;
            
            %Update the objects
            set(obj.(nam), 'BeamHandle', obj);
            set(obj.(nam), {'BeamEta'} , num2cell(objEta)');
            
        end        
        function removeBeamObject(obj, beamObj, eta)
            %removeBeamObject Removes a beam property object from a 
            %particular non-dimensional position along the beam.
            
            error('Update code');
            
        end
        function replaceBeamObject(obj, beamObj, eta)
            %replaceBeamObject Replaces the current distribution of beam
            %property objects of type 'class(beamObj)' with the beam
            %objects 'beamObj' at points 'eta' along the span.
            
            %Let 'assignBeamObject' handle the details...
            assignBeamObject(obj, beamObj, eta, 'replace');
                        
        end      
        % retrieving a specific beam property
        function BP = getBeamProperty(obj, tag)
            %getBeamProperty Return the handle to the specific beam
            %property objects with Quantity matching 'name'.
            
            %Screen for direct match with 'Quantity'
            allQuantities = {obj.BeamProperties.Quantity};  
            idx           = ismember(allQuantities, tag);
            
            %Maybe the user provided the 'Type' instead?
            if nnz(idx) == 0
               allTypes = {obj.BeamProperties.Type}; 
               idx      = ismember(allTypes, tag);
            end
            
            %Grab the beam properties 
            BP  = obj.BeamProperties(idx);
            
            if isempty(BP) %Force empty matrix
                BP = [];
            end
            
        end        
        % shortcut functions for 'awi.mixin.BeamProperty'
        function val = getEta(obj, varargin)
            %getEta Shortcut for calculating the finest possible
            %distribution of points along the notional beam axis.
            
            if isempty(obj.BeamProperties) %Check for empty
                val = [];
                return
            end
            
            %Pass it on
            val = getEta(obj.BeamProperties, varargin{:});
        end
        function val = getBPV(obj, Quantity, varargin)
            %getBPV Shortcut for calculating the beam property 'Quantity'
            %at a specified distribution along the notional beam.
            
            if nargin < 2 %Quit if no beam property provided
                val = [];
                return
            end
            
            if isempty(obj.BeamProperties) %Check for empty
                val = [];
                return
            end
            
            %Pass it on
            val = getBPV(obj.BeamProperties, Quantity, varargin{:});
        end        
    end
            
    methods (Access = private, Sealed) %Callbacks for dynamic property listeners
        %Default 'BeamProperty' behaviour
        function updateBeamGeometry(obj, src, ~)
            %updateBeamGeometry Passes the current value of 'XData',
            %'YData' or 'ZData' to the underlying beam properties so that
            %the 'BeamProperty' objects have an accurate understanding of
            %the beam geometry.
            
            %What has been changed and what is the new value?
            prpName = src.Name;
            prpVal  = obj.(prpName);
            
            %Propogate the changes to the underlying beam properties
            switch prpName
                case 'XData'
                    set(obj.BeamProperties, 'x', prpVal)
                case 'YData'
                    set(obj.BeamProperties, 'y', prpVal)
                case 'ZData'
                    set(obj.BeamProperties, 'z', prpVal)
                otherwise
                    error(['Unknown beam geometry has been passed to ', ...
                        'the listener method ''updateBeamGeometry''. ', ...
                        'Check the constructor of the ''Beamable'' '  , ...
                        'class.']);
            end  
            
        end
        function updateBeamPropertyAxis(obj, src, ~)
            %updateBeamGeometry Passes the current value of 'AxisFlag' to 
            %the underlying beam properties that are of type 'Planform
            %Property'. This is only relevant to objects of type
            %'awi.model.LiftingSurface'.
            
            %What has been changed and what is the new value?
            prpName = src.Name;
            prpVal  = obj.(prpName);
            
            %Down-select to get only those properties that are unique to
            %the 'LiftingSurface'.
            BP = selectBeamProp(obj.BeamProperties, obj.BeamPropType);
            
            %Propogate the changes to the underlying beam properties
            set(BP, 'AxisFlag', prpVal);
            
            %Make the changes to all dynamic properties ('_flag') so that
            %the public facing properties are a true reflection of the
            %underlying beam properties
            prpNames  = strcat({BP.Quantity}, '_flag');
            prpValues = repmat({prpVal}, size(prpNames));
            set(obj, prpNames, prpValues);

        end
        function updateBeamProperty(obj, src, ~)
            %updateBeamProperty Propogates the changes in a dynamic beam
            %property into the relevant 'BeamProperty' object in the object
            %array 'obj.BeamProperties'.
            
            %Which dynamic property has been changed?
            dPrpName = src.Name;
            
            %Which beam property are we dealing with?
            [bpName, BeamProp] = retrieveBeamPropName(obj, dPrpName);
            
            %Find the new value
            newVal = obj.(dPrpName);
            
            %Find this dynamic property's equivalent 'BeamProperty' object
            BPNames = {obj.BeamProperties.Quantity};
            idx     = ismember(BPNames, bpName);
            
            %Update the value of the 'Value' property in the 'BeamProperty'
            obj.BeamProperties(idx).(BeamProp) = newVal;
            
            %Propogate any changes to the dynamic prop
            %   - This is in case the 'BeamProperty' object does any
            %     further processing/validation of the data
            obj.(dPrpName) = obj.BeamProperties(idx).(BeamProp);
            
        end
        function val = retrieveBeamProperty(obj, src, ~)
            
            %Which dynamic property has been changed?
            dPrpName = src.Name;
            
            %Retrieve beam property name
            ind    = strfind(dPrpName, '_');
            if numel(ind) > 1
                ind = ind(end);
            end
            bpName = dPrpName(1 : ind - 1);
            
            %Retrieve the interpolated value
            val = getBPV(obj, bpName);
            
            %Replace NaN with 0
            val(isnan(val)) = 0;
            
            %Assign it back to the object
            obj.(dPrpName) = val;

        end
        %Validating BeamPropertyObjects
        function validateBeamPropObj(obj, src, ~)
            
            %Which beam property object are we dealing with?
            cls    = obj.BeamObjectClass(ismember(obj.BeamObjectNames, src.Name));            
            newVal = obj.(src.Name);
            
            %Allow characters/empty for now so that we can assign UUID 
            %during import process.
            if ~ischar(newVal) && ~isempty(newVal)
               validateattributes(obj.(src.Name), cls, {'row'}, class(obj), src.Name); 
            end
        end
        %Spawned 'BeamProperty' behaviour        
        function deriveBeamProperty(obj, src, ~)
            %deriveBeamProperty Derives the beam property specified by
            %'src.Name' from the internally stored beam property objects.
            
            %Which property?
            pName = src.Name;
            
            %Determine which beam property object the property belongs to
            bpoName = retriveBeamPropObjName(obj, pName);
            
            if isempty(obj.(bpoName)) %Escape route
                return
            end
            
            %If we get this far then the beam property object is populated
            %and we can derive the requested property from the stored beam
            %property object.
            newVal      = [obj.(bpoName).(pName)];
            if isscalar(newVal)
                newVal = repmat(newVal, size(obj.([pName, '_eta'])));
            end
            obj.(pName) = newVal;
            
            %Update the internally stored 'BeamProperty' object
            BP = getBeamProperty(obj, pName);
            BP.Value = newVal;
            
        end
        function updateSpawnedBeamProperty(obj, src, ~)
            %updateSpawnedBeamProperty Updates the 'Distribution',
            %'Variation' & 'AxisFlag' property values for the spawned beam
            %properties related to the beam proprty object which has
            %'src.Name' as one of its dynamic properties.
            
            %Which dynamic property has been changed?
            dPrpName = src.Name;
            
            %What is the value?
            newVal = obj.(dPrpName);
            
            %Which beam property object does this relate to?
            [bpoName, ~, sfx] = retrieveBeamPropName(obj, dPrpName);
            
            %What spawned dynamic beam properties belong to this beam
            %property object?
            idx = ismember(obj.BeamObjectNames, bpoName);
            bpNames = obj.BeamObjectProps{idx};
            
            %Assign the value to the spawned beam properties
            bpNames = cellfun(@(x) [x, sfx], bpNames, 'Unif', false);
            newVal  = repmat({newVal}, size(bpNames));
            set(obj, bpNames, newVal);
            
        end
    end 
   
    methods (Access = private, Sealed) %Helper functions
        function p = makeDynProps(obj, name, val)
            %makeDynamicProps Returns the handle to the dynamic properties
            %that have been added to the object with the property name
            %generated from 'name' and the suffixes defined by
            %'obj.DynPropTag'. The value is set as 'val'.
            
            %Create the dynamic properties which will depend on the beam
            %property object
            p = arrayfun(@(i)addDynamic(obj, [name, obj.DynPropTag{i}], ...
                val{i}), 1 : numel(val), 'Unif', false);
            p = horzcat(p{:});
            
        end
        function [bpName, BeamProp, sfx] = retrieveBeamPropName(obj, dPrpName)
            %retrieveBeamPropName Returns the name of the beam property and
            %the name of the underlying property in the
            %'awi.model.BeamProperty' object for a given dynamic property 
            %name.
            %- TODO : Update this to "contains" once we are sure of ver.
                        
            %Use logical indexing to allow for expandable list of beam 
            %properties            
            idx = cellfun(@(x) numel(x) >0 && ...
                strncmp(fliplr(dPrpName), fliplr(x), numel(x)), obj.DynPropTag);
            
            %Down select
            if any(idx)
                %Trim at the underscore ('_')
                BeamProp = obj.DynPropKey{idx};
                ind      = strfind(dPrpName, obj.DynPropTag{idx});
                bpName   = dPrpName(1 : ind - 1);
                sfx      = dPrpName(ind : end);
            else                
                BeamProp = 'Value';
                bpName   = dPrpName;
                sfx      = '';
            end
            
        end
        function bpoName = retriveBeamPropObjName(obj, pName)
            %retriveBeamPropObjName Returns the beam property object that
            %has the property with name 'pName'.
            
            %Check the stored list of beam property objects and their
            %relevant properties.
            idx = cellfun(@(x) any(ismember(x, pName)), obj.BeamObjectProps);
            assert(any(idx), sprintf(['No beam property object with ', ...
                'name ''%s''.'], pName));
            assert(nnz(idx) == 1, sprintf(['Expected only one beam ' , ...
                'property object with name ''%s''.'], pName));
            bpoName = obj.BeamObjectNames{idx};
            
        end
        function [nam, eNam] = parseBeamPropInputs(obj, beamObj, eta)
            %parseBeamPropInputs Parses the inputs to the beam property
            %object method and returns the name of the property(ies)
            %corresponding to this beam property object type.
                    
            %Sensible defaults
            nam  = [];
            eNam = [];
            
            %Only valid for scalar objects
            assert(numel(obj) == 1, ['Function ''assignBeamObject'' ', ...
                'not valid for handle arrays.']);
            
            %The beam property objects should all be of the same class and
            %at the very least should be of class 'awi.model.Entity'.
            assert(isa(beamObj, 'awi.model.Entity'), ['Expected the ', ...
                'beam property object to be of class ''awi.model.Entity.']);
            cls = unique(arrayfun(@(i) class(i), beamObj, 'Unif', false));
            assert(numel(cls) == 1, ['Expected the beam ', ...
                'property objects to be of the same class']);
            
            %The eta distribution should be a row vector of values
            %between [0 : 1]
            validateattributes(eta, {'numeric'}, {'row', '>=', 0, '<=', ...
                1}, class(obj), 'assignBeamObject');    
            
            %Check the object is a valid beam property object
            idx  = ismember(obj(1).BeamObjectClass, cls);
            if ~any(idx)
                warning(['Unable to assign a beam property object of ', ...
                    'class ''%s'' to the object as this is not a '    , ...
                    'valid beam property object for ''%s''.'],  cls{1} , ...
                    obj.NameAndType);
                return
            end
            
            %Which beam property object does this relate to?
            nam  = obj(1).BeamObjectNames{idx};
            eNam = [nam, '_eta'];
            
        end
        
    end
    
end

