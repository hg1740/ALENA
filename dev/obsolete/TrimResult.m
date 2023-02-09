classdef (ConstructOnLoad) TrimResult < awi.model.ResultSet
    %
    % Represents the result of performing a Trim analysis on a Aircraft

    properties % (SetAccess = immutable)
        
        %Extend metadata with content specific to awi
        LoadCase;
        RotMat;
        ControlDeflection;
        Loads;
        
    end
    
    properties (Dependent)
        
        Fx
        Fy
        Fz
        Mx
        My
        Mz
        Eta
        
        Legend
        
    end
    
    methods %get / set
        
        function val = get.Legend(obj)
           
            val = {obj.Name};
            
        end
        
        function val = get.Eta(obj)
            
           
            val = linspace(0,1,numel(obj.Loads{1}(1:6:end)))';
            
        end
        
        function val = get.Fx(obj)
            
            val = obj.Loads{1}(1:6:end);
            
        end
        
        function val = get.Fy(obj)
            
            val = obj.Loads{1}(2:6:end);
            
        end
        
        function val = get.Fz(obj)
            
            val = obj.Loads{1}(3:6:end);
            
        end
        
        function val = get.Mx(obj)
            
            val = obj.Loads{1}(4:6:end);
            
        end
        
        function val = get.My(obj)
            
            val = obj.Loads{1}(5:6:end);
            
        end
        
        function val = get.Mz(obj)
            
            val = obj.Loads{1}(6:6:end);
            
        end
        
    end
    
    methods % construction / destruction
        
        function obj = TrimResult(varargin)
            %
            % Construct a TrimResult object OBJ, based on the Trim Analysis PAR
        
            %Pass it on
            obj@awi.model.ResultSet(varargin{:});
            
            %Extend property groups
            obj.addPropertyGroup('General', ...
                'LoadCase', 'LoadCase');
            
        end
        
    end
    
end
