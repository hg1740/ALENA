classdef DrawableThing < mvc.model.Thing & mvc.mixin.Drawable
    %
    % Extend Thing class with functionality from Drawable
    
    methods % construction / destruction
        
        function obj = DrawableThing(varargin)
        
            %Pass it on
            obj@mvc.model.Thing(varargin{:});
            
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
                fcn = @mvc.model.DrawableThing;
    
            end
            
            %Start with parent
            obj = fcn(varargin{:});
           
            %For each child
            for i = 1:n
                
                %Add child
                ch = obj.add(fcn);
            
                %Give it a position, so it appears in different place from others when drawn
                ch.Position = [1, 0, 0] .* i;
                
                %Extend child
                ch.add(mvc.model.DrawableThing.defaultList(varargin{:}));
            
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
                fcn = @mvc.model.DrawableThing;
    
            end
            
            %Loop on size
            for i = 1:sz(1)
                for j = 1:sz(2)
                    
                    %Create instance
                    obj(i,j) = fcn(varargin{:});
                    
                    %Give it a position, so it appears in different place from others when drawn
                    obj(i,j).Position = [0, i, j];
                    
                end
            end
            
        end
        
    end
    
end
