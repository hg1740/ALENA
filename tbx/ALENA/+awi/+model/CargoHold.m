classdef CargoHold < awi.model.Compartment
    %CargoHold Defines a compartment that contains cargo. 
    %
    % Basically a thin wrapper around the 'awi.model.Compartment' object. 
    
    properties
        %Density of the cargo/payload [kgm^3]
        PayloadDensity
    end
    
    methods % set / get
        function set.PayloadDensity(obj, val)  %set.PayloadDensity 
            %set.PayloadDensity Set method for the property 
            %'PayloadDensity'.
            %
            % Passes the value straight to the 'Compartment' mass
            % density.
            
            obj.MassDensity = val;
            
        end
        function val = get.PayloadDensity(obj) %get.PayloadDensity 
           %get.PayloadDensity Get method for the property 
           %'PayloadDensity'.
           %
           % Retrieves the 'Compartment' mass density.
           
           val = obj.MassDensity;
           
        end
    end
    
    methods % construction
        function obj = CargoHold
           %CargoHold Constructor for the 'awi.model.CargoHold' object.
           %
           % Actions performed:
           %    - Sets the value of the 'CompartmentColour' and
           %    'PayloadColour' properties.
           
           obj.CompartmentColour = [18 , 35 , 224] ./ 255;
           obj.PayloadColour     = [68 , 252, 243] ./ 255;
           
        end
    end
    
end

