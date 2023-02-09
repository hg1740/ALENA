classdef FuelTank < awi.model.Compartment
    %FuelTank Defines a compartment that contains fuel. 
    %
    % Basically a thin wrapper around the 'awi.model.Compartment' object. 
    
    properties
        %Density of the fuel [kgm^3]
        FuelDensity
    end
    
    methods % set / get
        function set.FuelDensity(obj, val)  %set.FuelDensity 
            %set.FuelDensity Set method for the property 'FuelDensity'.
            %
            % Passes the value straight to the 'Compartment' mass
            % density.
            
            obj.MassDensity = val;
            
        end
        function val = get.FuelDensity(obj) %get.FuelDensity 
           %get.FuelDensity Get method for the property 'FuelDensity'.
           %
           % Retrieves the 'Compartment' mass density.
           
           val = obj.MassDensity;
           
        end
    end
    
    methods % construction
        function obj = FuelTank
           %FuelTank Constructor for the 'awi.model.FuelTank' object.
           %
           % Actions performed:
           %    - Sets the value of the 'CompartmentColour' and
           %    'PayloadColour' properties.
           
           obj.CompartmentColour = [227, 178, 43] ./ 255;
           obj.PayloadColour     = [161, 22 , 10] ./ 255;
           
        end        
    end
    
end

