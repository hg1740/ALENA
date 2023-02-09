classdef (ConstructOnLoad) FuelMass < awi.model.PointMass
    %FuelMass Represents a lumped mass of fuel.
    %
    % This class is a very thin wrapper around the 'awi.model.PointMass'
    % object and it only really exists in order for fuel mass to be grouped
    % together in the tree view.
    
    methods % constructor
        function obj = FuelMass(varargin)
            
            %Pass on to superclass
            obj@awi.model.PointMass(varargin{:});
            
            %Plot the point mass using a different marker style.
            obj.MarkerStyle     = 'd';
            obj.MarkerFaceColor = 'r';
            obj.MarkerEdgeColor = 'k';
            
        end
    end
    
end

