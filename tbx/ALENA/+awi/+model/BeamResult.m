classdef (ConstructOnLoad) BeamResult < awi.model.ResultSet & awi.mixin.Beamable
    %BeamResult Defines a sets of results over a beam.
    %
    % The results quantities are defined at specified non-dimensional
    % positions along a chosen axis of the beam. For more information see
    % the documentation for 'awi.mixin.Beamable'.
    %
    % An object of class 'awi.model.BeamResult' has the following results
    % quantites:
    %
    %   Internal Forces:
    %   - Fx : Internal load in the local beam x-axis (N)
    %   - Fy : Internal load in the local beam y-axis (N)
    %   - Fz : Internal load in the local beam z-axis (N)
    %   - Mx : Internal moment about the local beam x-axis (Nm)
    %   - My : Internal moment about the local beam y-axis (Nm)
    %   - Mz : Internal moment about the local beam z-axis (Nm)
    %
    %   Displacements:
    %   - Tx : Translation in the global x-direction
    %   - Ty : Translation in the global y-direction
    %   - Tz : Translation in the global z-direction
    %   - Rx : Rotation in the global x-direction
    %   - Ry : Rotation in the global y-direction
    %   - Rz : Rotation in the global z-direction
    %
    % These results quantities are added as dynamic properties using the
    % 'awi.mixin.Beamable' and 'awi.mixin.BeamProperty' functionality.
    %
    %   - TODO : Add rotation matricies as a dynamic property
        
    %Handle to parent beam
    properties 
        %Handle to the 'awi.model.Beam' object that these results are
        %related to
        Beam        
    end
    
    properties
        %Rotation Matrix describing the deformed to global translation at
        %certain eta positions.
        RotationMatrix 
    end
    
    %Useful shortcuts
    properties (Dependent)
        %Name of the beam that these results relate to
        BeamName
    end
    
    methods % set / get
        function set.Beam(obj, val)      %set.Beam     
            if isempty(val)
                obj.Beam = [];
                return
            end
            validateattributes(val, {'awi.model.Beam'}, {'scalar'}, ...
                class(obj), 'Beam');
            obj.Beam = val;
        end
        function val = get.BeamName(obj) %get.BeamName 
            %get.BeamName Get method for the dependent property 'BeamName'.
            
            %Pretty simple...
            if isempty(obj.Beam)
                val = '';
            else
                val = obj.Beam.Name;
            end
            
        end
    end
    
    methods % construction
        
        function obj = BeamResult(varargin)
            %BeamResult Constructor for the class 'awi.model.BeamResult'.
            
            %Default args
            args = { ... 
                'IsLeafNode'        , true, ... Result Sets inherrently have no children (? always true ?)
                'StructureLocked'   , true, ... follows implicitly from IsLeafNode, not really necessary
                'PropertiesLocked'  , true, ... do not allow user to edit content
                'PropertiesLockable', false ... user is not granted direct access to change editability
                };
            
            %Pass it on
            obj@awi.model.ResultSet(args{:}, varargin{:});
                        
            %Add results quantities...
            
            %Internal Loads
            addBeamProperty(obj, 'Fx', 'Type', 'Internal Loads', ...
                'Units', 'N');
            addBeamProperty(obj, 'Fy', 'Type', 'Internal Loads', ...
                'Units', 'N'); 
            addBeamProperty(obj, 'Fz', 'Type', 'Internal Loads', ...
                'Units', 'N'); 
            addBeamProperty(obj, 'Mx', 'Type', 'Internal Loads', ...
                'Units', 'Nm'); 
            addBeamProperty(obj, 'My', 'Type', 'Internal Loads', ...
                'Units', 'Nm'); 
            addBeamProperty(obj, 'Mz', 'Type', 'Internal Loads', ...
                'Units', 'Nm'); 
            
            %Internal loads in the FAME coordinate system
            addBeamProperty(obj, 'FxFAME', 'Type', 'Internal Loads', ...
                'Units', 'N', 'Name', 'Fx (FAME)');
            addBeamProperty(obj, 'FyFAME', 'Type', 'Internal Loads', ...
                'Units', 'N', 'Name', 'Fy (FAME)');
            addBeamProperty(obj, 'FzFAME', 'Type', 'Internal Loads', ...
                'Units', 'N', 'Name', 'Fz (FAME)');
            addBeamProperty(obj, 'MxFAME', 'Type', 'Internal Loads', ...
                'Units', 'Nm', 'Name', 'Mx (FAME)');
            addBeamProperty(obj, 'MyFAME', 'Type', 'Internal Loads', ...
                'Units', 'Nm', 'Name', 'My (FAME)');
            addBeamProperty(obj, 'MzFAME', 'Type', 'Internal Loads', ...
                'Units', 'Nm', 'Name', 'Mz (FAME)');
            
            %Displacements
            addBeamProperty(obj, 'Tx', 'Type', 'Displacements', ...
                'Name', 'Global-X Displacement', 'Units', 'm');
            addBeamProperty(obj, 'Ty', 'Type', 'Displacements', ...
                'Name', 'Global-Y Displacement', 'Units', 'm');
            addBeamProperty(obj, 'Tz', 'Type', 'Displacements', ...
                'Name', 'Global-Z Displacement', 'Units', 'm');
            addBeamProperty(obj, 'Rx', 'Type', 'Displacements', ...
                'Name', 'Global-X Rotation', 'Units', 'm');
            addBeamProperty(obj, 'Ry', 'Type', 'Displacements', ...
                'Name', 'Global-Y Rotation', 'Units', 'm');
            addBeamProperty(obj, 'Rz', 'Type', 'Displacements', ...
                'Name', 'Global-Z Rotation', 'Units', 'm');
            
            %Stresses?
            
            
        end
        
    end
    
end

