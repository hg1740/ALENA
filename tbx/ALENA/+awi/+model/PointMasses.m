classdef (ConstructOnLoad) PointMasses < awi.model.Collector
    %PointMasses Collector for a group of 'awi.model.PointMass' objects.
    %
    % Handles the plotting of 'awi.model.PointMass' objects in a vectorised
    % fashion.
    %
    % Responds to the 'buildElement' method and updates the name to reflect
    % the total mass of the child objects.
    %
    % Parameter Sets
    %   - 'mSet' (Mass Set) : Thin parameter set to provide a way to
    %   programatically change the name of the Collector to include the
    %   total mass of the collection.
    %
        
    properties (Dependent)
        %Type of mass that is being collected
        MassType
    end
    
    methods % set / get
        function val = get.MassType(obj) %get.MassType
            %get.MassType Get method for the dependent property 'MassType'.

            val = 'Point Masses';
            
            if obj.HasChildren && ...
                    isa(obj.Children(1), 'awi.model.FuelMass')
                val = 'Fuel Masses';
            end
            
        end
    end
    
    methods % construction
        function obj = PointMasses(varargin)
            
            %Pass on to superclass
            obj@awi.model.Collector(varargin{:});
            
            %Add some parameter sets
            obj.addParameterSet('mSet', ...
                'DisplayName', 'Mass Set', ...
                'Description', ['Allows the name of the collector to be ', ...
                'updated to reflect how much mass it contains'], ...
                'Precedence', inf, ...
                'MassType'  , 'Type of mass being collected.');
            
        end
    end
    
    methods % visualisation
        function hg = drawElement(obj, ht, tag)
            
            if nargin < 3
                tag = obj.MassType;
            end
            
            %Grab all point masses collected by this object
            PointMasses = obj.Children;
            
            if isempty(PointMasses) %Escape route
                hg{1} = gobjects(1);
                return
            end
            
            %Get the position of all point masses that are collected by
            %this object
            pos = get(PointMasses, 'Position');
            
            if iscell(pos) %Force matrix format for vectorised plotting                
                pos = vertcat(pos{:});
            end            
                        
            %Properties that might be specified
            prp = {'MarkerStyle', 'MarkerSize', 'MarkerFaceColor', 'MarkerEdgeColor'};
            val = get(PointMasses(1), prp);
            
            %Allow for propname mismatch
            prp{1} = 'Marker';
            
            %Throw away anything not specified
            prp(cellfun(@isempty, val)) = [];
            val(cellfun(@isempty, val)) = [];
            
            %Combine
            pav = [prp; val];
            
            %Draw the point masses 
            hg{1} = line('Parent', ht, ...
                'XData'    , pos(:, 1), ...
                'YData'    , pos(:, 2), ...
                'ZData'    , pos(:, 3), ...
                'Tag'      , tag      , ...
                'LineStyle', 'none'   , ... %Force discrete markers
                pav{:});
        end
    end
    
    methods % class building
        function build_mSet(obj)
            %build_mSet Updates the collector name to reflect the amount of
            %mass in the collection.
            
            %If we have any children
            if obj.HasChildren
                
                %Total mass of children - Round up to nearest whole number   
                M = ceil(sum([obj.Children.Mass]));
        
            else
                
                %No mass
                M = 0;
                
            end
            
            %Update the name to include masses 
            obj.Name = [obj.MassType, sprintf(' (%ikg)', M)];
        end
    end
    
end

