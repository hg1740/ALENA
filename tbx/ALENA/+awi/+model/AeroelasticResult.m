classdef (ConstructOnLoad) AeroelasticResult < awi.model.BeamResult
    %AeroelasticResult Defines a set of aeroelastic results over a beam.
    %
    % 'awi.model.AeroelasticResult' is of class 'awi.model.BeamResult',
    % therefore, it contains all results quantities in the superclass as
    % well as the following additional quantities:
    %
    %   - FxAero : External aerodynamic load in the ...
    %   - FyAero :
    %   - FzAero :
    %   - MxAero :
    %   - MyAero :
    %   - MzAero :
    %
    %   - Cl     : Lift coefficient
    %   - Cl_l   : Local lift coefficient (TODO - Update the description)
    %   - Cm     : Pitching moment coefficient
    
    properties
        ControlSurfaceDeflections
    end
    
    methods % construction
        function obj = AeroelasticResult(varargin) 
           
            %Pass it on
            obj@awi.model.BeamResult(varargin{:});
            
            %Add results quantities...
            
            %Aerodynamic forces in the global coordinate system
            %   - Lift is in the positive z-direction
            %   - Drag is in the positive x-direction (streamwise)
            addBeamProperty(obj, 'Lift', 'Type', 'Aerodynamic Loads');
            
            %Aeroydnamic forces in the beam local coordinate system
            addBeamProperty(obj, 'FxAero', 'Type', 'Aerodynamic Loads', ...
                'Name', 'Fx (Aero)', 'Label', 'F_{x_AERO}', ...
                'Description', 'The aerodynamic force in the beam local x-coordinate system');
            addBeamProperty(obj, 'FyAero', 'Type', 'Aerodynamic Loads', ...
                'Name', 'Fy (Aero)', 'Label', 'F_{y_AERO}', ...
                'Description', 'The aerodynamic force in the beam local y-coordinate system'); 
            addBeamProperty(obj, 'FzAero', 'Type', 'Aerodynamic Loads', ...
                'Name', 'Fz (Aero)', 'Label', 'F_{z_AERO}', ...
                'Description', 'The aerodynamic force in the beam local z-coordinate system'); 
            addBeamProperty(obj, 'MxAero', 'Type', 'Aerodynamic Loads', ...
                'Name', 'Mx (Aero)', 'Label', 'M_{x_AERO}', ...
                'Description', 'The aerodynamic moment about the beam local x-coordinate system'); 
            addBeamProperty(obj, 'MyAero', 'Type', 'Aerodynamic Loads', ...
                'Name', 'My (Aero)', 'Label', 'M_{y_AERO}', ...
                'Description', 'The aerodynamic moment about the beam local y-coordinate system'); 
            addBeamProperty(obj, 'MzAero', 'Type', 'Aerodynamic Loads', ...
                'Name', 'Mz (Aero)', 'Label', 'M_{z_AERO}', ...
                'Description', 'The aerodynamic moment about the beam local z-coordinate system'); 
            
            %Aerodynamic Force & Moment Coefficients
            addBeamProperty(obj, 'Cl'    , 'Type', 'Aerodynamic Loads', ...
                'Name', 'Lift Coefficient', 'Label', 'C_{L}', ...
                'Description', 'The lift coefficient');
            addBeamProperty(obj, 'Cl_l'  , 'Type', 'Aerodynamic Loads', ...
                'Name', 'Normalised Lift Coefficient (Cl_l)', 'Label', 'C_{L_l}');
            addBeamProperty(obj, 'Cm'    , 'Type', 'Aerodynamic Loads', ...
                'Name', 'Pitching Moment', 'Label', 'C_{M}', ...
                'Description', 'The pitching coefficient');
            
        end
    end
    
end

