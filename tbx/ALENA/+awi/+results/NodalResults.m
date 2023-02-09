classdef NodalResults < matlab.mixin.SetGet
    %NodalResults Defines a set of translation and rotation results
    %quantities 
    
    %Results data
    properties
        %Handle to the 'awi.fe.Node' objects -> nNode = numel(obj.Nodes)
        Nodes  
        %Time vector - > nT = numel(obj.TimeVector)
        TimeVector      
        %Translational terms in global frame [3, nNode, nT]
        Translation 
        %Rotational terms [3, nNode, nT]
        Rotation        
    end
    
    %Descriptors
    properties
        %Results type
        ResultsType
        %Title describing the results SET
        Title
    end
    
    properties (Constant)
        ValidResultsType = {'Displacement', 'Velocity', 'Acceleration', 'Eigenvector'};
    end
    
    methods % set / get
        function set.Nodes(obj, val)        %set.Nodes
            validateattributes(val, {'awi.fe.Node', 'awi.fe.NodeCollection'}, ...
                {'row'}, class(obj), 'Nodes');
            obj.Nodes = val;
        end
        function set.Translation(obj, val)  %set.Translation
            validateattributes(val, {'numeric'}, {'3d', 'nrows', 3, ...
                'finite'}, class(obj), 'Translation');
            obj.Translation = val;
        end
        function set.Rotation(obj, val)     %set.Rotation
            validateattributes(val, {'numeric'}, {'3d', 'nrows', 3, ...
                'finite'}, class(obj), 'Rotation');
            obj.Rotation = val;
        end
        function set.TimeVector(obj, val)   %set.TimeVector
            validateattributes(val, {'numeric'}, {'row', 'nonnan', ...
                'finite', 'real', 'nonnegative'}, class(obj), 'TimeVector');
            obj.TimeVector = val;
        end
        function set.ResultsType(obj, val)  %set.ResultsTypes
            val = validatestring(val, obj.ValidResultsType, class(obj), 'ValidResultsTypes');
            obj.ResultsType = val;
        end
        function set.Title(obj, val)        %set.Title
            validateattributes(val, {'char'}, {'row'}, class(obj), 'Title');
            obj.Title = val;
        end
    end
    
end

