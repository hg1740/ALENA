classdef Properties < mvc.view.Container

    properties
        
        %The model to which we are attached
        Model;
        
        %The field to display
        Field = 'Name';
        
        %Function to call on selection change
        SelectionChangeFcn;
        
        %Track currently selected item in model
        % (required if model is a collection, because this view can present only one item)
        Selection;
        
        %Control whether we respond to externally applied changes in selection
        SelectionLinked = true;
        
    end
    
    properties (Access = protected)
        
        %Selection popup
        hSelection;
        
        %SelectionLinked checkbox
        hSelectionLinked;
        
        %The panel in which we display properties
        hPanel;
        
        %Placeholder for listeners
        Listeners;  
        
    end
    
    methods % get/set
        
        function val = get.Selection(obj)
            
            %Get from uicontrol
            val = get(obj.hSelection, 'Value');
            
            %Return the actual selection, rather than the index of the selection
            val = obj.Model(val);
            
        end
        
        function set.Selection(obj, val)
        
            %Interested in responding to externally applied change in selection ?
            if ~obj.SelectionLinked %#ok<MCSUP>
                
                %No
                return;
                
            end
               
            %Make a note
            obj.Selection = val;
            
            %No selection is valid from user's point of view
            if isempty(val)
                
                %But we need something to put into popup
                idx = 1;
                
            elseif isnumeric(val)
                
                %TODO: ensure within range ?
                idx = val;
                
            elseif isa(val, class(obj.Model)) %#ok<MCSUP>

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
                mdl = obj.Model; %#ok<MCSUP>
                b = arrayfun(@(x)handle(val) == handle(mdl(x)), 1:numel(mdl));
                
                %TODO: this test not quite adequate ?
                assert(any(b), 'invalid selection');
                
                %Make a note
                idx = find(b);
                
            end
            
            %Assign in uicontrol
            set(obj.hSelection, 'Value', idx); %#ok<MCSUP>
            
            %Force an update
            update(obj);
            
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
            
            %Make a note
            obj.Model = val;
            
            %Update content
            update(obj);
            
        end
        
    end
    
    methods % construction
        
        function obj = Properties( model, varargin )
            
            %To allow the constructor to be called with no arguments (e.g. for test purposes)
            if nargin == 0
                
                %Make something up
                model = mvc.model.DrawableThing.defaultTree;
                
            end
            
            %Call superclass constructor
            obj@mvc.view.Container(varargin{:});
            
            %Start with a VBox
            hv = uiextras.VBox('Parent', obj.UIContainer);
            
            %Containing a grid
            hg = uiextras.Grid('Parent', hv);
            
            %Containing a label, and popup (to control selection)
            uicontrol('Parent', hg, 'Style', 'Text', 'String', 'Selection:')
            obj.hSelection = uicontrol('Parent', hg, 'Style', 'popup', 'callback', @obj.onSelectionChange);
            hg.Widths = [75, -1];
            
            %Add a checkbox for the 'SelectionLinked' flag
            obj.hSelectionLinked = uicontrol('Parent', hg, ...
                'Style', 'Checkbox', ...
                'String', 'linked', ...
                'TooltipString', 'Selection in this View is linked to selection in Application', ...
                'Callback', @obj.onSelectionLinked);
            
            %And a panel
            %hSP = uix.ScrollingPanel('Parent', hv); %Edit 15/06/2018 - C.Szczyglowski
            obj.hPanel = uiextras.Panel('Parent', hv, 'Padding', 6);
            
            %Set heights
            hv.Heights = [25, -1];
            
            %Store the model
            obj.Model = model;
            
            %Listen for change
            obj.Listeners = event.listener(obj.Model, 'ModelChanged', @obj.onModelChanged);
            
        end
        
    end
    
    methods ( Access = private )
        
        function onModelChanged(obj, ~, ~ )
            
            %Update content
            update(obj);
            
        end
        
        function onSelectionChange(obj, ~, ~)
            
            %Update content
            update(obj);
                
            %Anything else to do ?
            if ~obj.SelectionLinked
                
                %No
                
            elseif ~isempty(obj.SelectionChangeFcn)
                
                %Pass it on
                obj.SelectionChangeFcn(obj.Selection);
                
            end
            
        end
        
        function onSelectionLinked(obj, ~, ~)
            
            %Make a note
            obj.SelectionLinked = boolean(get(obj.hSelectionLinked, 'Value'));
            
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
            sel = get(obj.hSelection, 'Value');
                
            %In case multiple selection has been specified (from external source)
            sel(2:end) = [];
            
            %Ensure selection popup remains valid
            sel = min(sel, numel(str));
            
            %In case of no selection
            if isempty(sel)
                set(obj.hSelection, 'String', {'no selection'}, 'Value', 1, 'Enable', 'off');
            else
                set(obj.hSelection, 'String', str, 'Value', sel, 'Enable', 'on');
            end
            
            %Linked option
            set(obj.hSelectionLinked, 'Value', obj.SelectionLinked);
            
            %Does selected item support context menus ?
            if isa(mdl(sel), 'mvc.mixin.Contextable')
                
                %Yes - pass it on
                context(mdl(sel), obj.hSelection);
                
            end            
            
            %Redraw view of selected element
            propedit(mdl(sel), obj.hPanel);
            
        end
        
    end
    
end
