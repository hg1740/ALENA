classdef (ConstructOnLoad) Shapeable < mvc.mixin.Drawable
    %
    % Represents a thing with shape, with a position and orientation in 3-D space
    %  inherited from Drawable
    
    properties (AbortSet, SetObservable)
        
        %Thing can be drawn as a shape in 3-D space
        Shape = [];
        
        %With dimensions set by the following
        Length = 1;
        Width = 1;
        Height = 1;
        
        %Drawn with this many facets
        Facets = 20;
        
        %Tagged with this
        Tag = 'Shape';
        
        %Appearance of the shape
        FaceColor;
        FaceAlpha;
        Edges = false;
        
    end
    
    properties (Dependent)
        
        Radius;
        
    end
    
    properties (SetAccess = protected, Transient)
        
        Shapes = {'sphere', 'hemisphere', 'cylinder', 'hemicylinder', 'cuboid'};
        
    end
    
    methods % get/set
    
        function set.Radius(obj, val)
        
            %What shape are we ?
            switch lower(obj.Shape)
                
                case {'sphere', 'hemisphere'}
                    
                    %Apply to all
                    obj.Length = val;
                    obj.Width = val;
                    obj.Height = val;
                    
                case {'cylinder', 'hemicylinder'}
                    
                    %Apply to some
                    
                otherwise
                    
            end
                                
        end
        
        function val = get.Shape(obj)
            
            %Start here
            val = obj.Shape;
            
            %Nothing yhet ?
            if isempty(val)
                
                %Use this by default
                val = obj.Shapes{1};
                
            end
            
        end
        
    end
    
    methods % construction / destruction
        
        function obj = Shapeable(varargin)
            
            %Start with base class
            obj@mvc.mixin.Drawable(varargin{:});
            
            %If properties are managed by Nameable
            if isa(obj, 'mvc.mixin.Nameable')
                
                %Extend property groups
                obj.addPropertyGroup('Geometry', ...
                    'Shape', 'Shape', [], obj.Shapes, ...
                    'Length', 'Length', [], [], ...
                    'Width', 'Width', [], [], ...
                    'Height', 'Height', [], [], ...
                    'Tag', 'Tag', [], [], ...
                    'Facets', 'Facets', [], []);

                obj.addPropertyGroup('Appearance', ...
                    'FaceColor', 'Face color', [], [], ...
                    'FaceAlpha', 'Face alpha', [], [], ...
                    'Edges', 'Edges', [], []);
                
                %Geometry and Appearance only relevant if object is Visible
                obj.setPropertyGroup('Geometry', 'Visible', @obj.Visible);
                obj.setPropertyGroup('Appearance', 'Visible', @obj.Visible);
                
            end
            
        end
        
    end
    
    methods % visualisation
        
        function hg = drawElement(obj, ht, tag)
            
            %Caller supplied tag explicitly ?
            if nargin < 3
                
                %No - apply default tag
                tag = obj.Tag;
                
            end
            
            %Start with base class
            hg = drawElement@mvc.mixin.Drawable(obj, ht);
            
            %Get x,y,z data carefully (e.g. so user can edit content bit by bit)
            [xd, yd, zd] = xyzdata(obj);
            
            %Get additional arguments
            args = {'FaceColor', obj.FaceColor; ...
                'FaceAlpha', obj.FaceAlpha; ...
                'EdgeColor', mvc.mixin.UiTools.bool2offon(obj.Edges, {'none', []}); ...
                'Tag', tag};
            
            %Eliminate any empties
            args(cellfun(@isempty, args(:,2)),:) = [];
            
            %Transpose before apply
            args = args.';
            
            %Draw the shape
            hg{end+1} = surf('Parent', ht, ...
                'XData', xd, ...
                'YData', yd, ...
                'ZData', zd, ...
                args{:});
            
        end
        
    end
    
    methods (Access = private)
        
        function [xd, yd, zd] = xyzdata(obj)
            
            %What are we ?
            switch lower(obj.Shape)
                
                case {'sphere', 'hemisphere'}
                    
                    %Pass it on
                    [xd, yd, zd] = sphere(obj.Facets);
                    
                    %Scale up
                    xd = xd .* obj.Length;
                    yd = yd .* obj.Width;
                    zd = zd .* obj.Height;
                    
                    %For hemisphere
                    if strcmpi(obj.Shape, 'hemisphere')
                        
                        %Throw half away
                        n = ceil(size(xd,2)/2) + 1;
                        xd(:,n:end) = [];
                        yd(:,n:end) = [];
                        zd(:,n:end) = [];
                        
                    end
                    
                case {'cylinder', 'hemicylinder'}
                    
                    %Allow for stretch
                    n = unique([numel(obj.Length), numel(obj.Width)]);
                    assert(numel(unique(n(n > 1))) <= 1, 'dimension mismatch');

                    %Pass it on
                    [xd, yd, zd] = cylinder(ones(1, max(n)), obj.Facets);
                    
                    %Scale up
                    xd = xd .* obj.Length(:);
                    yd = yd .* obj.Width(:);
                    zd = zd .* obj.Height(:);
                    
                    %For hemicylinder
                    if strcmpi(obj.Shape, 'hemicylinder')
                        
                        %Throw half away
                        n = ceil(size(xd,2)/2) + 1;
                        xd(:,n:end) = [];
                        yd(:,n:end) = [];
                        zd(:,n:end) = [];
                        
                    end
                 
                case 'cuboid'
                    
                    %Pass it on - a cube is same as cylinder with 4 facets
                    [xd, yd, zd] = cylinder(1, 4);
                    
                    %Scale up, the factor of sqrt(2) required to true-up the dims
                    % (otherwise the unit length is on the diagonal)
                    xd = xd .* obj.Length .* sqrt(2);
                    yd = yd .* obj.Width .* sqrt(2);
                    zd = zd .* obj.Height .* sqrt(2);
                    
                    %Followed by a rotation
                    M = makehgtform('zrotate',deg2rad(45));
                    
                    %Strip the last row and col
                    M(end,:) = [];
                    M(:,end) = [];
                    
                    %Splurge
                    xyz = [xd(:), yd(:), zd(:)];
                    
                    %Rotate
                    xyz = xyz * M;
                    
                    %Unsplurge
                    xd = reshape(xyz(:,1), size(xd));
                    yd = reshape(xyz(:,2), size(yd));
                    zd = reshape(xyz(:,3), size(zd));                    
                    
                otherwise
                    error(['shape ''', obj.Shape, ''' not (yet) handled']);
            end
            
        end
        
    end
    
end
