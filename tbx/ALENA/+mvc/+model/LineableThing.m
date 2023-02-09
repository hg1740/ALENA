classdef LineableThing < mvc.model.Thing & mvc.mixin.Lineable
    %
    % Extend Thing class with functionality from Lineable
    
    methods % construction / destruction
        
        function obj = LineableThing(varargin)
        
            %Pass it on
            obj@mvc.model.Thing(varargin{:});
            
        end
        
    end
    
    methods (Static)
        
        function obj = simpleExample(varargin)
            
            %Constructor can be specified explicitly
            if nargin > 1 && isa(varargin{1}, 'function_handle')
                
                %Get from inputs
                fcn = varargin{1};
                varargin(1) = [];
                
            else
                
                %Just go with this
                fcn = @mvc.model.LineableThing;
    
            end
            
            %Start with parent
            obj = fcn(varargin{:});
           
            %With line data as follows
            obj.XData = [0, 100];
            
            %Add content
            ch1 = obj.add(fcn);
            ch1.Offset = 0.4;
            ch1.XData = [0, 50];
            ch1.Rotation = [0, 0, 60];
            
            %Add content
            ch2 = obj.add(fcn);
            ch2.Offset = 0.4;
            ch2.XData = [0, 50];
            ch2.Rotation = [0, 0, -60];
            
            %Add content
            ch3 = obj.add(fcn);
            ch3.Offset = 1;
            ch3.XData = [0, 10];
            ch3.Rotation = [0, -75, 0];
            
            %Add content
            ch31 = ch3.add(fcn);
            ch31.Offset = 1;
            ch31.XData = [0, 5];
            ch31.Rotation = [0, 30, 90];
            ch31.RotationOrder = [3, 2, 1];
            
            %Add content
            ch32 = ch3.add(fcn);
            ch32.Offset = 1;
            ch32.XData = [0, 5];
            ch32.Rotation = [0, 30, -90];
            ch32.RotationOrder = [3, 2, 1];
            
        end
        
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
                fcn = @mvc.model.LineableThing;
    
            end
            
            %Start with parent
            obj = fcn(varargin{:});
           
            %For each child
            for i = 1:n
                
                %Add child
                ch = obj.add(fcn);
            
                %Give it some data
                ch.XData = [0, 1] .* i;
                ch.YData = [0, 0];
                ch.ZData = [0, 0];
                
                %Give it a position, so it appears in different place from others when drawn
                ch.Position = [1, 0, 0] .* i;
                
                %Extend child
                ch.add(mvc.model.LineableThing.defaultList(varargin{:}));
            
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
                fcn = @mvc.model.LineableThing;
    
            end
            
            %Loop on size
            for i = 1:sz(1)
                for j = 1:sz(2)
                    
                    %Create instance
                    obj(i,j) = fcn(varargin{:});
            
                    %Give it some data
                    obj(i,j).XData = [0, 1] .* i;
                    obj(i,j).YData = [0, 0];
                    obj(i,j).ZData = [0, 0];
                    
                    %Give it a position, so it appears in different place from others when drawn
                    obj(i,j).Position = [0, i, j];
                    
                end
            end
            
        end
        
    end
    
end
