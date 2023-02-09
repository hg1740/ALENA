classdef (ConstructOnLoad) Collector < mvc.model.Collector & awi.mixin.Buildable
    %Collector Represents a node in a hierarchical model that serves no 
    %purpose other than to collect child nodes and support parametric class
    %definitions.
    %
    % On construction the Collector suppresses any detail of property 
    % groups and removes propedit from context.
    %
    % The 'awi.model.Collector' class is subtlely different to the
    % 'mvc.model.Collector' class in that it provides the functionality for
    % collectors to implement methods belonging to the
    % 'awi.mixin.Buildable' class. 
    %
    % This means AWI-Collectors can respond to the 'ModelChanging' event by
    % implementing 'Parameter Sets' (see 'awi.mixin.Buildable') which will
    % allow objects to be updated on the fly!
    %
    % TODO : The object hierachy is a bit messy when we inherit from
    % 'mvc.model.Collector' as now there are two ways into the 'mvc'
    % package. Ideally I would prefer a single point of entry (i.e.
    % awi.model.Entity') but as the Collector classes have specific
    % functions I think this is okay...
    
    methods % construction
        function obj = Collector(varargin)
            
            %Start with base class
            obj@mvc.model.Collector(varargin{:});
            
            %Reinstate the property groups to maintain compatibility with
            %'mvc.mixin.Drawable' (See constructor of 'mvc.mixin.Drawable'
            %for more information)
            obj.addPropertyGroup('Geometry');
            obj.setPropertyGroup('Geometry', 'Visible', @obj.Visible);
            obj.addPropertyGroup('General' , 'Visible', 'Visible');
            
            %Set the collector as visible so that the 'drawElement' method
            %is triggered
            obj.Visible = true;
            
        end
    end
    
end

