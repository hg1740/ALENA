classdef (ConstructOnLoad) PointMass < awi.model.Component
    %PointMass Defines a point mass with inertia properties.
    
    %Mass & Inertia properties
    properties (AbortSet, SetObservable)
        %Value of the mass acting a point
        Mass      = 0;
        %Inertia in the 11 direction as specified by ...
        Inertia11 = 0;
        Inertia22 = 0;
        Inertia33 = 0;
        Inertia12 = 0;
        Inertia23 = 0;
        Inertia13 = 0;
    end
    
    %Mass Group 
    properties
        %Defines a sub-set of the model mass (e.g. Primary, Secondary,
        %Systems, etc.) This will be used to tag the masses in the
        %LayeredDrawing view.
        MassGroup = '';
    end
    
    %Inertia matrix
    properties (Dependent)
        InertiaMatrix        
    end
    
    methods % set / get
        function set.Mass(obj, val)           %set.Mass          
            %set.Mass Set method for the property 'Mass'.
            %
            % 'Mass' must be a scalar numeric
            
            validateattributes(val, {'numeric'}, {'scalar', 'nonempty', ...
                'finite', 'real', 'nonnan'}, class(obj), 'Mass');
            obj.Mass = val;
        end
        function set.Inertia11(obj, val)      %set.Inertia11     
            %set.Inertia11 Set method for the property 'Inertia11'.
            %
            % 'Inertia11' must be a scalar numeric
            
            validateattributes(val, {'numeric'}, {'scalar', 'nonempty', ...
                'finite', 'real', 'nonnan'}, class(obj), 'Inertia11');
            obj.Inertia11 = val;
        end
        function set.Inertia22(obj, val)      %set.Inertia22     
            %set.Inertia22 Set method for the property 'Inertia22'.
            %
            % 'Inertia22' must be a scalar numeric
            
            validateattributes(val, {'numeric'}, {'scalar', 'nonempty', ...
                'finite', 'real', 'nonnan'}, class(obj), 'Inertia22');
            obj.Inertia22 = val;
        end
        function set.Inertia33(obj, val)      %set.Inertia33     
            %set.Inertia33 Set method for the property 'Inertia33'.
            %
            % 'Inertia33' must be a scalar numeric
            
            validateattributes(val, {'numeric'}, {'scalar', 'nonempty', ...
                'finite', 'real', 'nonnan'}, class(obj), 'Inertia33');
            obj.Inertia33 = val;
        end
        function set.Inertia12(obj, val)      %set.Inertia12     
            %set.Inertia12 Set method for the property 'Inertia12'.
            %
            % 'Inertia12' must be a scalar numeric
            
            validateattributes(val, {'numeric'}, {'scalar', 'nonempty', ...
                'finite', 'real', 'nonnan'}, class(obj), 'Inertia12');
            obj.Inertia12 = val;
        end
        function set.Inertia23(obj, val)      %set.Inertia23     
            %set.Inertia23 Set method for the property 'Inertia23'.
            %
            % 'Inertia23' must be a scalar numeric
            
            validateattributes(val, {'numeric'}, {'scalar', 'nonempty', ...
                'finite', 'real', 'nonnan'}, class(obj), 'Inertia23');
            obj.Inertia23 = val;
        end   
        function set.Inertia13(obj, val)      %set.Inertia13     
            %set.Inertia12 Set method for the property 'Inertia13'.
            %
            % 'Inertia13' must be a scalar numeric
            
            validateattributes(val, {'numeric'}, {'scalar', 'nonempty', ...
                'finite', 'real', 'nonnan'}, class(obj), 'Inertia13');
            obj.Inertia13 = val;
        end
        function set.MassGroup(obj, val)      %set.MassGroup     
            %set.MassGroup Set method for the property 'MassGroup'
            %
            %   - 'MassGroup' must be a character row vector.
            
            if isempty(val) %Force character
                obj.MassGroup = '';
                return
            end
            validateattributes(val, {'char'}, {'row', 'vector', 'nonempty'}, ...
                class(obj), 'MassGroup');            
            obj.MassGroup = val;
        end
        function val = get.InertiaMatrix(obj) %get.InertiaMatrix 
            %get.InertiaMatrix Get method for the dependent property
            %'InertiaMatrix'.
            %
            % 'InertiaMatrix' is the 3x3 matrix of mass moment of inertia
            % values for this point mass. 
            %
            % It is assumed that the inertia matrix is symmetric.
            
            val = [ ...
                obj.Inertia11, obj.Inertia12, obj.Inertia13 ; ...
                obj.Inertia12, obj.Inertia22, obj.Inertia23 ; ...
                obj.Inertia13, obj.Inertia23, obj.Inertia33];
            
        end
    end
    
    methods % constructor
        function obj = PointMass(varargin)
           
            %Pass on to superclass
            obj@awi.model.Component(varargin{:});
            
            %Update property groups
            obj.addPropertyGroup('Mass and Inertia', ...
                'Mass', 'Point Mass', 'Mass acting at a point', [], ...
                'Inertia11', 'I11'  , 'Inertia in the plane-11', []);
            
            %Plot the point mass using a different marker style.
            obj.MarkerStyle     = 's';
            obj.MarkerFaceColor = 'b';
            obj.MarkerEdgeColor = 'k';
            
            %PointMass has no children
            obj.IsLeafNode = true;
            
            %Use a custom collector object
            obj.CollectorClass = 'awi.model.PointMasses';
            
            %TODO - Finesse this!
            obj.Visible = false;
                        
        end
    end
    
end

