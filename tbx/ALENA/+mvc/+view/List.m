classdef List < mvc.view.Container

    properties
        
        %The model to which we are attached
        Model;
        
        %The field to display in list
        Field = 'Name';
        
        %Function to call on selection change
        SelectionChangeFcn;
        
        %Support multiple selection ?
        AllowMultipleSelect = false;
        
    end
    
    properties (Dependent)
        
        %Track currently selected item in model
        Selection;
        
    end
    
    properties (Access = protected)
        
        %The listbox in which we display the model
        hSelection;
        
        %Placeholder for listeners
        Listeners;
        
    end
    
    methods
        
        function val = get.Selection(obj)
            
            %Get from uicontrol
            val = get(obj.hSelection, 'Value');
            
            %Ensure selection is valid
            n = numel(obj.Model);
            val(val < 1 | val > n) = [];
            
            %Nothing left ?
            if isempty(val)
                
                %Is that a valid response ?
                if ~obj.AllowMultipleSelect
                    
                    %No
                    val = 1;
                    
                end
                
            end
            
            %Return the actual selection, rather than the index of the selection
            val = obj.Model(val);
            
        end
        
        function set.Selection(obj, val)
            
            %Cater for no selection
            if isempty(val)
                
                %Nothing to select in listbox, if multiple-select is OK
                if obj.AllowMultipleSelect
                    idx = [];
                else
                    
                    %Need to show something
                    idx = 1;
                    
                end
            
            elseif isnumeric(val)
                
                %TODO: assert in range
                idx = val;
                
            elseif isa(val, class(obj.Model))

                %Allow selection to be set with an actual member of the list
                %
                %As long as it is in the list
                % Would like to do this, but it fails
                %[b, idx] = ismember(val, obj.Model);
                %assert(all(b), 'invalid selection');
                %
                %This doesn't work either, under all circumstances, because of heterogeneous collection
                %idx = arrayfun(@(x)find(x == obj.Model), val);
                %assert(all(idx > 0), 'invalid selection');
                
                %Instead
                mdl = obj.Model; 
                b = arrayfun(@(x)handle(val) == handle(mdl(x)), 1:numel(mdl));
                
                %Make a note
                idx = find(b);
                
            end
            
            %Assign in uicontrol
            set(obj.hSelection, 'Value', idx);
            
        end
        
        function val = get.Model(obj)
            
            %Get the underlying value
            val = obj.Model;
            
            %If Collectable
            if isa(val, 'mvc.mixin.Collectable')
                
                %Flatten it
                val = flatlist(val, 2);
                
            end
            
        end
        
        function set.Model(obj, val)
            
            %Allow model to be nothing at al
            if isequal(val, [])
                
            else
                
                %Valid class ?
                assert(isa(val, 'mvc.mixin.Nameable') || isa(val, 'mvc.mixin.Collectable'), 'model MUST be Collectable or Nameable');

            end
            
            %Make a note
            obj.Model = val;
            
            %Update content
            update(obj);

        end
        
        function obj = List(model, varargin )
            
            %To allow the constructor to be called with no arguments (e.g. for test purposes)
            if nargin == 0
                
                %Make something up
                model = mvc.model.DrawableThing.defaultTree;
                
            end
            
            %Call superclass constructor
            obj@mvc.view.Container( varargin{:} );
            
            %Initialise listbox
            obj.hSelection = uicontrol('Parent', obj.UIContainer, ...
                'Units', 'normalized', ...
                'Position', [0, 0, 1, 1], ...
                'Style', 'listbox', ...
                'Min', 0, ...
                'Max', 1 + obj.AllowMultipleSelect, ...
                'Callback', @obj.onSelectionChange);
            
            %Store the model
            obj.Model = model;
            
            %Listen for change
            obj.Listeners = event.listener(obj.Model, 'ModelChanged', @obj.onModelChanged);
            
        end
        
    end
    
    methods (Access = private)
        
        function onModelChanged(obj, ~, ~ )
            
            %Update content
            update(obj);
            
        end
        
        function onSelectionChange(obj, ~, ~)
            
            %Update content
            update(obj);
            
            %Anything else to do ?
            if ~isempty(obj.SelectionChangeFcn)
                
                %Pass it on
                obj.SelectionChangeFcn(obj.Selection);
                
            end
            
        end
        
        function update(obj)
            
            %Get the model
            mdl = obj.Model;
            
            %Anything to do ?
            if isempty(mdl)
                
                %Bail
                return;
                
            end
            
            %Get content to be displayed
            str = {mdl.(obj.Field)}; 
            
            %If the underlying list is actually hierarchical
            if isa(mdl, 'mvc.mixin.Collectable')
                
                %Then indent each member of string by depth of corresponding element in tree structure
                n = [mdl.DepthInTree];
                str = arrayfun(@(i)[repmat(' ', 1, n(i)), str{i}], 1:numel(str), 'UniformOutput', false);
                
            end
            
            %Get current selection
            sel = obj.Selection;
            
            %Does selected item support context menus ?
            if isa(sel, 'mvc.mixin.Contextable')
                
                %Yes - pass it on
                context(sel, obj.hSelection);
                
            end            
            
            %Verify membership
            [b, val] = ismember(sel, mdl);
            assert(b, 'selection not a member of list !');
            
            %Display in list
            set(obj.hSelection, 'String', str, 'Value', val);

        end
        
    end
    
end
