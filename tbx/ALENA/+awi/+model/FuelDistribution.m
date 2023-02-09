classdef (ConstructOnLoad) FuelDistribution < awi.model.Entity
    %FuelDistribution Handles a relationship between an 'awi.model.Beam'
    %object and a collection of 'awi.model.PointMass' objects.
    
    properties
        %Handle to the 'awi.model.Beam' object that these point masses are
        %associated with
        ParentBeam
        %Vector of 'awi.model.FuelMass' objects
        FuelMasses
    end
    
    properties (Hidden) %TODO - Get rid of these properties eventually
        %Vector of 'Conm2' objects from the AeroFlex module
        %   - TODO : This is just a temporary fix to get the fuel
        %   distribution working correctly. Get rid of eventually.
        AeroFlexConm2       
        %Name of the fuel file from FAME
        FameFuelFile
        %Character describing the aircraft part that these masses are
        %associated with
        PartName
    end
    
    methods % set / get
        function set.ParentBeam(obj, val) %set.ParentBeam
            %set.ParentBeam Set method for the dependent property
            %'ParentBeam'.
            %
            %   - 'ParentBeam' must be a scalar instance of the
            %     'awi.model.Beam' class.
            
            validateattributes(val, {'awi.model.Beam'}, {'scalar', ...
                'nonempty'}, class(obj), 'ParentBeam');
            obj.ParentBeam = val;
        end
        function set.FuelMasses(obj, val) %set.FuelMasses
            %set.FuelMasses Set method for the dependent property
            %'FuelMasses'.
            %
            %   - 'FuelMasses' must be a row vector of 'awi.model.FuelMass'
            %     objects.
            
            validateattributes(val, {'awi.model.FuelMass'}, {'row', ...
                'vector', 'nonempty'}, class(obj), 'FuelMasses');
            obj.FuelMasses = val;
        end
    end
    
    methods % construction
        function obj = FuelDistribution(varargin)
            %FuelDistribution Constructor of the class 'FuelDistribution'.
            %
            % TODO - Implement a generic method for pulling out inputs to
            % an object and assigning them to the object after calling the
            % constructor.
            
            %Caller may be supplying properties we need to grab at this level
            prp = {'ParentBeam', 'FuelMasses', 'AeroFlexConm2', 'FameFuelFile', 'PartName'};
            
            %Look for match in varargin
            idx = zeros(size(varargin));
            b = cellfun(@ischar, varargin);
            [~, idx(b)] = cellfun(@(x)ismember(x, prp), varargin(b));
            
            %Anything ?
            if any(idx)
                
                %Pull out corresponding values
                val = varargin(find(idx > 0) + 1);
                
                %Elliminate from varargin
                varargin(sort([find(idx > 0), find(idx > 0) + 1])) = [];
                
            end
            
            %So the properties we found in varargin are
            prp = prp(idx(idx > 0));
            
            %Pass it on to superclass
            obj@awi.model.Entity(varargin{:});
            
             %Assign values at this level, if any were passed in
            if ~isempty(prp)
                pav = [prp; val];
                set(obj, pav{:});
            end

        end
        
    end
    
    methods % attaching a fuel distribution to a AWI model
%         function attach
%             
%         end
    end
    
end

