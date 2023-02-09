classdef Thing < mvc.mixin.Nameable ... so each thing has metadata that are manageable
        & mvc.mixin.Searchable      ... so each thing, and any arrays of things, are searchable
        & mvc.mixin.Contextable     ... so each thing can present context-sensitive options to the user
        & mvc.mixin.Dynamicable     ... so each thing can have additional properties dynamically added if required
        & mvc.mixin.Collectable     ... so things can be arranged in a hierarchical fashion (with parent/child relationships)
        & mvc.mixin.UiTools         ... so things can interact with the user
        ... & mvc.mixin.Auditable       ... so things can maintain a log of changes, actions, etc
        & mvc.mixin.Debugable         % so the user can be presented with debug functions
    %
    % Base-class of all application specific objects comprising an array of generic things.
    
    methods % construction / destruction
        
        function obj = Thing(varargin)
        
            %Pass it on - ensuring contructors are called in correct order
            obj@mvc.mixin.Nameable(varargin{:});
            obj@mvc.mixin.Searchable;
            obj@mvc.mixin.Contextable;
            obj@mvc.mixin.Collectable;
            obj@mvc.mixin.UiTools;
            %obj@mvc.mixin.Auditable;
            obj@mvc.mixin.Debugable;

        end
        
    end
    
    methods (Static)
        
        function obj = defaultTree(n, varargin)
            
            %How many children ?
            if nargin < 1
                
                %Go with three
                n = 3;
                
            end
            
            %Constructor can be specified explicitly
            if nargin > 1 && isa(varargin{1}, 'function_handle')
                
                %Get from inputs
                fcn = varargin{1};
                varargin(1) = [];
                
            else
                
                %Just go with this
                fcn = @mvc.model.Thing;
    
            end
            
            %Start with parent
            obj = fcn(varargin{:});
           
            %For each child
            for i = 1:n
                
                %Add child
                ch = obj.add(fcn);
            
                %Extend child
                ch.add(mvc.model.Thing.defaultList(varargin{:}));
            
            end
            
        end
        
        function obj = defaultList(sz, varargin)
            
            %How many items in list ?
            if nargin < 1
                
                %Go with five
                sz = [5, 1];
                
            elseif numel(sz) == 1
                
                %Assume a column
                sz(2) = 1;
                
            end
            
            %Constructor can be specified explicitly
            if nargin > 1 && isa(varargin{1}, 'function_handle')
                
                %Get from inputs
                fcn = varargin{1};
                varargin(1) = [];
                
            else
                
                %Just go with this
                fcn = @mvc.model.Thing;
    
            end
            
            %Loop on size
            for i = 1:sz(1)
                for j = 1:sz(2)
                    
                    %Create instance
                    obj(i,j) = fcn(varargin{:});
                    
                end
            end
            
        end
        
    end
    
    methods (Access = protected) % to make dynamic properties copyable
        
        function cpy = copyElement(obj)
            
            %Start with dynamicable copy
            cpy = copyElement@mvc.mixin.Dynamicable(obj);
            
            %Then nameable copy
            cpy = copyElement@mvc.mixin.Nameable(cpy);
            
            %Then collectable copy
            cpy = copyElement@mvc.mixin.Collectable(cpy);
            
        end
        
    end
            
end
