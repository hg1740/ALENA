classdef (ConstructOnLoad) Buildable < handle
    %Buildable Allows an object to be defined parametrically using an
    %arbitrary number of parameterisation schemes.
    %
    % Implements a generic 'build' method that invokes 'buildElement' which
    % builds the object before calling 'build' on all children of the
    % object. This allows the hierachy of the model to be respected and by
    % ensuring that the parent is built before the child it is possible for
    % the child object to be dependent on the parent.
    %   e.g. consider the example of a strut-braced wing. To connect the
    %   strut to the wing the wing must first be built before the strut.
    %
    % Fundamental to the 'Buildable' class is the use of 'ParameterSets'. A
    % parameter set is a group of properties that can be used to completly
    % define all of the primary properties of the object. The 'Buildable'
    % class allows multiple parameter sets to be defined and then the user
    % can choose which parameter set should be used to defined the model.
    %
    % When parameter sets are added to the object using the
    % 'addParameterSet' method a property group is also defined containing
    % all of the properties in the parameter set.
    
    % Defining a parameter set
    % Title        - Must be of type char or string.
    % DisplayName  - Must be of type char or string.
    % Description  - Must be of type char or string
    % Precedence   - Must be a double which is a whole number
    % Parameters   - Must be a two column cell array of char or string
    % BuildMethod  - Must be a function handle
    % UpdateMethod - Must be a function handle
    properties (SetAccess = protected, Transient)
        ParameterSets = struct('Title', {} , 'DisplayName', {}, ...
            'Description' , {}        , 'Precedence' , {}, ...
            'Parameters'  , {{[], []}}, 'BuildMethod', {}, ...
            'UpdateMethod', {});
    end
    
    %ActiveSet, ActiveSetMode
    properties (AbortSet, SetObservable)
        %Defines which parameter set is used to build the model
        ActiveSet = 'none';        
        %Determines whether the model updates automatically if a parameter is changed
        ActiveSetMode = 'auto';        
    end
    
    properties (SetAccess = 'protected')        
        %A status message, internal maintained
        BuildStatus = '';        
    end
    
    %SetNames, ActiveSet_, AllParamNames, ActiveSetDisplayName
    properties (Dependent)
        %Cell array containing all the titles and display names of the parameter sets
        SetNames
        %Token used by the classdef to determine active set
        ActiveSet_
        %Cell array containing the name of all the parameters
        AllParamNames        
        %Helper
        ActiveSetDisplayName        
    end
    
    methods % set / get 
        function set.ActiveSetDisplayName(obj, val) %set.ActiveDisplayName
            
            %Pass it on
            obj.ActiveSet = val;
            
        end
        function set.ActiveSetMode(obj, val)        %set.ActiveSet
            %set.ActiveSetMode Set method for property 'ActiveSetMode'.
            %
            % Rules :
            %   - 'ActiveSetMode' must be a string/char array matching one
            %     of the following tokens: 'manual', 'auto'
            
            validatestring(val, {'manual', 'auto'}, class(obj), 'ActiveSetMode');
            obj.ActiveSetMode = val;
        end
        function set.ActiveSet(obj, val)            %set.ActiveSet
            %set.ActiveSet Set method for property 'ActiveSet'.
            %
            % Rules :
            %   - 'ActiveSet' must be a string/char array matching either
            %     the 'Title' or 'DisplayName' of one of the parameter sets
            %
            % 06/03/18,pjr: Need to allow 'none' else basic tests like
            %  export/import fail
            validatestring(val, [obj.SetNames(:)', {'none'}]);  %#ok<MCSUP>
            
            %Allow caller to specify by name or description
            [b, pdx] = ismember(val, obj.SetNames(2,:));
            if b
                val = obj.SetNames{1,pdx};
            end
            
            %Assign in object
            obj.ActiveSet = val;
            
        end
        function val = get.SetNames(obj)            %get.SetNames
            %get.SetNames Get method for the dependant property 'SetNames'
            
            val = [{obj.ParameterSets.Title}; {obj.ParameterSets.DisplayName}];
            
            if isempty(val)
                val = {'none'; 'none'};
            end
        end
        function val = get.ActiveSet_(obj)          %get.ActiveSet_
            %get.ActiveSet_ Get method for the dependent property 'ActiveSet_'
            %
            % The property 'ActiveSet' can be either the parameter set
            % Title or the parameter set DisplayName, however, for ease of
            % use the code will only use Title to identify the sets and
            % grab properties of the sets.
            %
            % This method return the parameter set title for a given value
            % of 'ActiveSet'.
            
            %Start with existing/default value
            val = obj.ActiveSet;
            
            %Titles and display names for all sets
            setNames = obj.SetNames;
            
            %Check for empty 'setNames'
            if isempty(setNames)
                return
            end
            
            %Has DisplayName been provided?
            idx = ismember(setNames(2, :), obj.ActiveSet);
            
            %Define 'ActiveSet_'
            if any(idx)
                val = setNames(1, idx);
            else
                val = obj.ActiveSet;
            end
            
        end
        function val = get.AllParamNames(obj)       %get.AllParamNames
            %get.AllParamNames Get method for the dependent property
            %'AllParamNames'
            %
            % 'AllParamNames' is a cell array containing all of the
            %  parameter names without any repetitions.
            
            %Extract just the parameter names from each parameter set
            c = cellfun(@(x) x(:, 1), {obj.ParameterSets.Parameters}, 'Unif', false);
            
            %Make column cell array and find unique parameter names
            val = unique(cat(1, c{:}));
            
        end
        function val = get.ActiveSetDisplayName(obj)%get.ActiveDisplayName
            
            %Start here
            val = obj.ActiveSet;
            
            %Look it up
            [b, pdx] = ismember(val, obj.SetNames(1,:));
            if b
                val = obj.SetNames{2,pdx};
            end
            
        end
    end
    
    methods % construction / destruction
        function obj = Buildable(varargin) %Buildable
            
            %Add in check for 'mvc.mixin.Nameable' and
            %'mvc.mixin.Collectable'
            
            %Might be useful
            if isa(obj, 'mvc.mixin.Contextable') && isa(obj, 'mvc.mixin.Debugable')
                obj.addContext('|Debug>Build...', 'build');
            end
            
            %Update Property Groups
            obj.addPropertyGroup('Build Set', ... enough to warrant its own tab page (?)
                'ActiveSetDisplayName', 'Active Set',       'Defines which parameter set is used to build the model', @()obj.SetNames(2,:), ...
                'ActiveSetMode',        'Active Set Mode',  'Determines whether the active set is automatically updated or fixed', {'auto', 'manual'}, ...
                'BuildStatus',          'Build status',     'Status of this object after last call to ''build''', []);
            
            %This property group only relevant at least one parameter set defined
            obj.setPropertyGroup('Build Set', ...
                'Visible', @()size(obj.SetNames, 2) > 1);
            
            %Root object
            if isa(obj, 'awi.model.Framework')
                
                %Listens to itself, to cascade calls to "build"
                addlistener(obj, 'ModelChanging', @obj.build);
                
            end
            
        end
    end
    
    methods (Sealed) % build methods
        function build(obj, ~, ~)                               %build           
            %build Invokes the 'buildElement' method for this object and
            %then calls the 'build' method for any children of this object.
            %
            % This way the entire model can be built by invoking the
            % 'build' method for the top level object in the model
            % hierachy.
            %
            % TODO - Enable object arrays
    
            assert(numel(obj) == 1, 'Can only build one object at a time');
             
            %Build this object
            buildElement(obj)
            
            %Build all children
            for iChild = 1 : numel(obj.Children) %#ok<*MCNPN>
                if isa(obj.Children(iChild), 'awi.mixin.Buildable')
                    build(obj.Children(iChild));
                end
            end
            
        end
        function buildElement(obj)                              %buildElement    
            %buildElement Builds the object using the appropriate
            %parameters as defined by the 'obj.ActiveSet'. If the required
            %parameters have not been defined then an appropriate method is
            %selected based on 'obj.ActiveSetMethod' and the precedence of
            %the parameter sets.
            
            %TODO - Enable object arrays for 'build' method
            assert(numel(obj) == 1, 'Can only build one object at a time');
            
            %Clear out any extant status message
            obj.BuildStatus = '';
            
            %If 'obj.ParameterSets' is empty then no parameter sets have
            %been defined for this object so no further action needed.
            if isempty(obj.ParameterSets)
                return
            end
            
            %Careful
            try
                
                %Get current state of the object
                %prps   = allProperties(obj);
                %prpVal = get(obj, prps); %--> Too expensive
                
                %Get names of all parameter sets
                setNames = obj.SetNames(1, :);
                
                % check for common parameters in 'obj.ParameterSets'
                %   - Assumption is that the common parameters are needed by
                %     every set otherwise they wouldn't be common!
                idx = strcmpi(setNames, 'Common');
                if any(idx)
                    cProps = obj.ParameterSets(idx).Parameters(:, 1);
                else
                    cProps = {};
                end
                
                % if 'obj.ActiveSet' has not been defined then use first
                % parameter set but make sure 'obj.ActiveSetMode' = 'auto'.
                if isempty(obj.ActiveSet) || strcmpi(obj.ActiveSet, 'none')
                    obj.ActiveSet     = obj.ParameterSets(1).Title;
                    obj.ActiveSetMode = 'auto';
                end
                
                f_getSetProps = @(set) [obj.ParameterSets(ismember(setNames, set)).Parameters(:, 1); cProps];
                
                % get property names and values for 'obj.ActiveSet' - include
                % any common properties as these may also be required
                activeSetProps = f_getSetProps(obj.ActiveSet);
                %activeSetVal   = prpVal(ismember(prps, activeSetProps));
                activeSetVal   = get(obj, activeSetProps);
                
                %Parent can be empty!
                idx = ~ismember(activeSetProps, 'Parent');
                activeSetVal   = activeSetVal(idx); 
                
                % try and build the object
                if ~any(cellfun('isempty', activeSetVal)) && ~isempty(activeSetVal)
                    
                    %Pause listeners
                    listeners(obj, 'pause');
                    
                    % get function handle to build method for 'obj.ActiveSet'
                    f = getMethodHandle(obj, 'build', obj.ActiveSet);
                    f(obj); % build the object using this method
                    % which sets require updating?
                    setNames(ismember(setNames, {'Common', obj.ActiveSet})) = [];                    
                    % update all other sets
                    for iSet = 1 : numel(setNames)
                        f = getMethodHandle(obj, 'update', setNames{iSet});
                        f(obj);
                    end
                    
                    %Run the update method to catch any
                    %other events that need to happen.
                    updateModel(obj);
                    
                    %Resume listeners
                    listeners(obj, 'resume');
                    
                    %What happened ?
                    obj.BuildStatus = ['built using active set ''', obj.ActiveSet, ''''];
                    
                    %Return to invoking function
                    return
                else
                    switch obj.ActiveSetMode
                        
                        case 'auto'
                            % make a decision about which set to use
                            
                            % remove the current active set and the 'Common' set
                            %    - Assumption is that the common set on its own
                            %      is never sufficient to define the object.
                            %    - Remove the current active set as we know it
                            %      is not fully populated.
                            reducedSet = setNames(~ismember(setNames, ...
                                {'Common', obj.ActiveSet}));
                            
                            % find a fully populated set and build the model
                            for iSet = 1 : numel(reducedSet)
                                % grab property names and values for this set
                                idx      = ismember({obj.ParameterSets.Title}, reducedSet{iSet});
                                setProps = [obj.ParameterSets(idx).Parameters(:, 1); cProps];
                                setVal   = get(obj, setProps);
                                %setVal   = prpVal(ismember(prps, setProps));
                                if ~any(cellfun('isempty', setVal))
                                    
                                    %Pause listeners
                                    listeners(obj, 'pause');
                                    
                                    % build the object
                                    f = getMethodHandle(obj, 'build', reducedSet{iSet});
                                    f(obj);
                                    % update 'obj.ActiveSet'
                                    obj.ActiveSet = reducedSet{iSet};
                                    % which sets need updating?
                                    setNames(ismember(setNames, ...
                                        {'Common', obj.ActiveSet})) = [];
                                    % run update methods for all other sets
                                    for uSet = 1 : numel(setNames)
                                        f = getMethodHandle(obj, 'update', setNames{uSet});
                                        f(obj);
                                    end
                                    
                                    %Run the update method to catch any
                                    %other events that need to happen.
                                    updateModel(obj);
                                    
                                    %What happened ?
                                    obj.BuildStatus = ['built using auto-selected active set ''', obj.ActiveSet, ''''];
                                    
                                    %Resume listeners
                                    listeners(obj, 'resume');
                                    
                                    % return to invoking function
                                    return
                                end
                            end
                            
                            % If we get this far then the model has not been built need to inform the user.
                            obj.BuildStatus = sprintf('unable to build using any of the available sets ''%s'', not enough information provided', ...
                                strjoin(setNames, ''', '''));
                            
                        case 'manual'
                            
                            % Unable to choose a different build method as the
                            % user has explicitly specified the active set.
                            % -> Print a message so they know that the object
                            %    cannot be built.
                            obj.BuildStatus = sprintf('unable to build using set ''%s'', not enough information provided', ...
                                obj.ActiveSet);
                            
                    end
                    
                end
                
            catch err
                
                %Make a note
                obj.BuildStatus = sprintf('FAILED to build using set ''%s'', with error ''%s''', ...
                    obj.ActiveSet, err.message);
                
            end
            
        end
        function f = getMethodHandle(obj, methodType, setName)  %getMethodHandle 
            %getMethodHandle Grabs the function handle for a method of type
            %'methodType' for the parameter set 'setName'.
            
            validatestring(methodType, {'build', 'update'});
            
            % grab current parameter set
            paramSet = obj.ParameterSets(ismember({obj.ParameterSets.Title}, setName));
            
            % check for a predefined method in the parameter set
            switch methodType
                case 'build'
                    f = str2func([methodType, '_', setName]);
                case 'update'
                    f = paramSet.UpdateMethod;
            end
            
            % if a method has not been specified for this parameter set
            % then construct the default method name.
            if isempty(f)
                f = str2func([methodType, '_', setName]);
            end
        end
        function addParameterSet(obj, tit, varargin)            %addParameterSet 
            
            %Optional key words in 'varargin'
            tokens = {'DisplayName', 'Description', 'Precedence', ...
                'BuildMethod', 'UpdateMethod', 'Precedence'};
            
            %Does this parameter group already exist?
            tdx = find(strcmp(tit, {obj.ParameterSets.Title}), 1);
            
            if isempty(tdx)
                
                %No - extend with new set
                tdx = numel(obj.ParameterSets) + 1;
                obj.ParameterSets(tdx).Title = tit;
                
            end
            
            %Caller may optionally provide 'DisplayName', 'Description',
            %'Precedence', 'BuildMethod' or 'UpdateMethod'
            %   - search varargin for these tokens and parse.
            if ~isempty(varargin)
                for iT = 1 : numel(tokens)
                    id = find(strcmp(varargin, tokens{iT}), 1);
                    if ~isempty(id)
                        obj.ParameterSets(tdx).(tokens{iT}) = varargin{id+1};
                        varargin([id, id+1]) = [];
                    end
                end
            end
            
            %How many properties are we adding ?
            n = numel(varargin) / 2;
            
            %How many parameters are already in this set?
            nP = size(obj.ParameterSets(tdx).Parameters, 1);
            
            %Assume that any remaining data are Parameter-Names &
            %Parameter-Description pairs.
            obj.ParameterSets(tdx).Parameters(nP+1:nP+n, 1) = varargin(1:2:end);
            obj.ParameterSets(tdx).Parameters(nP+1:nP+n, 2) = varargin(2:2:end);
            
            %If a precedence has not been defined then give it a high value
            if isempty(obj.ParameterSets(tdx).Precedence)
                obj.ParameterSets(tdx).Precedence = inf;
            end
            
            %Add new parameter set to obj.PropertyGroups
            tit = obj.ParameterSets(tdx).DisplayName;
            c = obj.ParameterSets(tdx).Parameters';
            obj.addPropertyGroup(tit, c{:});
            
            %If this is not the 'Common' property group
            if ~strcmpi(tit, 'Common Set')
                
                %It is only relevant if this parameter set is selected
                obj.setPropertyGroup(tit, 'Visible', @()strcmp(obj.ActiveSetDisplayName, tit));
                
            end
            
            %Update the order of the parameter sets based on the precedence
            [~, I] = sort([obj.ParameterSets.Precedence]);
            obj.ParameterSets = obj.ParameterSets(I);
            
        end
        function Set = getParameterSet(obj, tit)                %getParameterSet 
            %Retrieves a parameter set by referencing the 'Title' property
            %group or the 'DisplayName' property group.
            
            titles   = {obj.ParameterSets.Title};
            dispName = {obj.ParameterSets.DisplayName};
            
            ind = strcmpi(titles, tit);
            
            if ~any(ind)
                % perhaps 'tit' is the parameter set display name?
                ind = strcmpi(dispName, tit);
            end
            
            % return the Parameter Set or provide user with a message
            if ~any(ind)
                fprintf(['No Parameter Set exists with ''Title'' or ', ...
                    '''DisplayName'' : ''%s''.\n'], tit);
            else
                Set = obj.ParameterSets(ind);
            end
            
        end
    end
   
    methods % build methods
        function updateModel(~) %updateModel
            %updateModel Runs a final set of updates to the model to
            %capture any events/actions that need to take place after the
            %individual update methods have been executed.
            %
            % At this level of the class hierachy this function does not
            % actually do anything. Let the specifics be handled by the
            % subclasses.
        end
    end
    
end