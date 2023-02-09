classdef (ConstructOnLoad) Lineable < mvc.mixin.Drawable
    %
    % Represents a drawable thing, with a position and orientation in 3-D space
    %  inherited from Drawable
    %
    % This subclass extends Drawable by adding a list of points (XData, YData, ZData)
    %  and associated appearance properties (color, width etc)
    
    properties (AbortSet, SetObservable)
        
        %Thing can be drawn as a line in 3-D space
        XData = [];
        YData = [];
        ZData = [];
        
        %Appearance of the line
        LineColor;
        LineStyle;
        LineWidth;
        MarkerStyle;
        MarkerSize;
        MarkerFaceColor;
        MarkerEdgeColor;
         
    end
    
    properties (Dependent)
       RData 
    end
    
    methods % setters / getters
        function val = get.RData(obj)
            %get.RData Get method for the property 'RData'.
            %
            % 'RData' is a row vector containing the straight line distance
            % along the 'XData', 'YData', 'ZData' coordinates.
            
            [xd, yd, zd] = xyzdata(obj);            
            dXYZ = diff(xd).^2 + diff(yd).^2 + diff(zd).^2;
            val  = cumsum([0, sqrt(dXYZ)]);
        
            if isscalar(val)
               val = [0, 1]; 
            end
        end
    end
    
    methods % construction / destruction
        
        function obj = Lineable(varargin)
            
            %Start with base class
            obj@mvc.mixin.Drawable(varargin{:});
            
            %If properties are managed by Nameable
            if isa(obj, 'mvc.mixin.Nameable')
                
                %Extend property groups
                obj.addPropertyGroup('Geometry', ...
                    'XData', 'x-data', ...
                    'YData', 'y-data', ...
                    'ZData', 'z-data');
                
                %Extend property groups
                obj.addPropertyGroup('Appearance', ...
                    'LineColor', 'Line color', ...
                    'LineStyle', 'Line style', ...
                    'LineWidth', 'Line width', ...
                    'MarkerStyle', 'Marker style', ...
                    'MarkerSize', 'Marker size', ...
                    'MarkerFaceColor', 'Marker face color', ...
                    'MarkerEdgeColor', 'Marker edge color');
            
                %Appearance only relevant if object is Visible
                obj.setPropertyGroup('Appearance', 'Visible', @obj.Visible);
                
                %Exporting appearance to file is only relevant if object is visible
                obj.setPropertyGroup('Appearance', 'Export', @obj.Visible);
                
            end
            
        end
        
    end
    
    methods % visualisation
        
        function hg = drawElement(obj, ht, tag)

            %Caller supplied tag explicitly ?
            if nargin < 3
                
                %No - apply a tag indicating that we are just drawing some lines
                tag = 'Lines';
                
            end
            
            %Start with base class
            hg = drawElement@mvc.mixin.Drawable(obj, ht);
            
            %Get x,y,z data carefully (e.g. so user can edit content bit by bit)
            [xd, yd, zd] = xyzdata(obj);
            
            %Get additional arguments
            args = {'Color', obj.LineColor; ...
                'LineStyle', obj.LineStyle; ...
                'LineWidth', obj.LineWidth; ...
                'Marker', obj.MarkerStyle; ...
                'MarkerSize', obj.MarkerSize; ...
                'MarkerFaceColor', obj.MarkerFaceColor; ...
                'MarkerEdgeColor', obj.MarkerEdgeColor; ...
                'Tag', tag};
            
            %Eliminate any empties
            args(cellfun(@isempty, args(:,2)),:) = [];
            
            %Draw the line
            hg{end+1} = line('Parent', ht, ...
                'XData', xd, ...
                'YData', yd, ...
                'ZData', zd, ...
                args{:});                
            
        end
        
        function pos = s2pos(obj, s)
            %
            %Given s in range [0, 1] return the coordinates [x, y, z]
            % of the point the fraction s along the length of the line
                
            %Start from x,y,z data
            [xd, yd, zd] = xyzdata(obj);
            
            %If we have nothing at all
            if isempty(xd)
                
                %Then send back all zeros (? or error ?)
                pos = [0, 0, 0];
                
            elseif s == 0
                
                %We're at the start
                pos = [xd(1), yd(1), zd(1)];
                
            elseif s == 1
                
                %We're at the end
                pos = [xd(end), yd(end), zd(end)];
                
            else
                
                %What's the length of each line segment ?
                ld = sqrt(sum([xd; yd; zd] .^ 2, 1));
                
                %So the cumulative length of the line ?
                lc = cumsum(ld);
                
                %How far along the line are we looking ?
                x = s .* lc(end);
                
                %Which line segment is this ?
                n = find(x > lc, 1, 'first');
                
                %What fraction along the line-segment ?
                f = 1 - (lc(n + 1) - x) ./ ld(n + 1);
                
                %So the point of interest is
                pos(1) = xd(n) + f .* diff(xd(n:n+1));
                pos(2) = yd(n) + f .* diff(yd(n:n+1));
                pos(3) = zd(n) + f .* diff(zd(n:n+1));
                
            end
        
        end
        
    end
    
    methods (Access = private)
    
        function [xd, yd, zd] = xyzdata(obj)
        
            %In principle it's as simple as this
            xd = obj.XData;
            yd = obj.YData;
            zd = obj.ZData;
            
            %BUT, to allow user to build up content incrementally, we can
            % ensure that the call to 'line' won't fall over just because
            % the content is work in progress
            n = cellfun(@numel, {xd, yd, zd});
            
            %Nothing yet ?
            if all(n == 0)
                
                %Nothing more to worry about
                return;
                
            end
            
            %Ensure we send back no more than the shortest non-empty list
            nmax = min(n(n > 0));
            xd(nmax+1:end) = [];
            yd(nmax+1:end) = [];
            zd(nmax+1:end) = [];
            
            %And if anything is unspecified, we can just go with zeros
            if isempty(xd)
                xd = zeros(1, nmax);
            end
            if isempty(yd)
                yd = zeros(1, nmax);
            end
            if isempty(zd)
                zd = zeros(1, nmax);
            end
            
        end
        
    end
    
end
