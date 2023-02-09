classdef AuditTrail < mvc.view.Container

    properties
        
        %The model to which we are attached
        Model;
        
        %Function to call on selection change
        SelectionChangeFcn;
        
        %Track currently selected item(s)
        Selection;
        
    end
    
    properties (Access = protected)
        
        %The table in which we display the audit trail
        hTable;
        
        %Placeholder for listeners
        Listeners;  
        
    end
    
    methods
        
        function set.Model(obj, val)
            
            %Make a note
            obj.Model = val;
            
            %Update content
            update(obj);
            
        end
        
        function obj = AuditTrail( model, varargin )
            
            %To allow the constructor to be called with no arguments (e.g. for test purposes)
            if nargin == 0
                
                %Make something up
                model = mvc.model.DrawableThing.defaultTree;
                
            end
            
            %Call superclass constructor
            obj@mvc.view.Container( varargin{:} );
            
            %Initialise table
            hPanel = uiextras.Panel('Parent', obj.UIContainer, 'Padding', 6);
            obj.hTable = uitable('Parent', hPanel, 'Units', 'normalized', 'Position', [0, 0, 1, 1]);            
            
            %Store the model
            obj.Model = model;
            
            %Listen for change
            obj.Listeners = event.listener(obj.Model, 'ModelChanged', @obj.onModelChanged);
            
        end
        
    end
    
    methods ( Access = private )
        
        function update(obj)
            
            %Get the model
            mdl = obj.Model;
            
            %If Collectable
            if isa(mdl, 'mvc.mixin.Collectable')
                
                %Flatten it
                mdl = flatlist(mdl, 2);
                
            end
            
            %Only interested in Auditable content
            b = arrayfun(@(x)isa(x, 'mvc.mixin.Auditable'), mdl);
            mdl = mdl(b);
            
            %Anything to do ?
            if isempty(mdl)
                
                %Bail
                return;
                
            elseif ~isa(mdl, 'mvc.mixin.Auditable')
                
                %Bail
                return;
                
            end
            
            %Get ALL audit trail details
            dat = vertcat(mdl.AuditTrail);
            
            %Cater for time-tables
            if (~verLessThan('matlab', '9.3') && istimetable(dat))
                
                %Sort the table (TODO: should be a better way than this)
                [t, tdx] = sort(dat.Properties.RowTimes, 'descend');
                dat = dat(tdx,:);
                
            else
                
                %Sort the table
                [t, tdx] = sort(dat.When, 'descend');
                dat = dat(tdx,:);
                
            end
           
            %If a time-table (need version check, for backward compat)
            if ~verLessThan('MATLAB', '9.3') && istimetable(dat)
                
                %Set row names in table
                set(obj.hTable, 'RowName', char(t));
                
            end
            
            %Display in table
            hdg = dat.Properties.VariableNames;
            dat = mvc.mixin.Nameable.tabulatable(table2cell(dat));
            set(obj.hTable, 'Data', dat, 'ColumnName', hdg);
            
        end
        
        function onModelChanged(obj, ~, ~ )
            
            %Update content
            update(obj);
            
        end
        
    end
    
end
