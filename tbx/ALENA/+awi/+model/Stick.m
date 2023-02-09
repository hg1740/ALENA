classdef (ConstructOnLoad) Stick < awi.model.Component
    %Stick Describes a line in 3D space.
    
    %Line data
    properties (AbortSet, SetObservable)
        %Thing can be drawn as a line in 3-D space
        XData = [];
        YData = [];
        ZData = [];
    end
       
    properties (Dependent) %RData, Length 
        %Straight line distance along the Stick
        RData
        %Maximum straight line distance along the 'Stick'
%         Length
    end
    
    %Appearance of the line
    properties
        LineColor;
        LineStyle;
        LineWidth;
    end
    
    methods % set / get
        function set.XData(obj, val)   %set.XData
            %set.XData Set method for the property 'XData'
            %
            %   - 'XData' must be a row vector with finite, real values.
            
            validateattributes(val, {'numeric'}, {'row', 'finite', ...
                'real'}, class(obj), 'XData');
            obj.XData = val;
        end
        function set.YData(obj, val)   %set.YData
            %set.YData Set method for the property 'YData'
            %
            %   - 'YData' must be a row vector with finite, real values.
            
            validateattributes(val, {'numeric'}, {'row', 'finite', ...
                'real'}, class(obj), 'YData');
            obj.YData = val;
        end
        function set.ZData(obj, val)   %set.ZData
            %set.ZData Set method for the property 'ZData'
            %
            %   - 'ZData' must be a row vector with finite, real values.
            
            validateattributes(val, {'numeric'}, {'row', 'finite', ...
                'real'}, class(obj), 'ZData');
            obj.ZData = val;
        end
        function val = get.RData(obj)  %get.RData
            %get.RData Get method for the property 'RData'.
            %
            % 'RData' is a row vector containing the straight line distance
            % along the 'XData', 'YData', 'ZData' coordinates.
            
            [xd, yd, zd] = xyzdata(obj);            
            val = obj.getLineLength(xd, yd, zd);
            if iscolumn(val)
                val = val';
            end
            if isscalar(val)
                val = [0, 1];
            end
        end
%         function val = get.Length(obj) %get.Length
%             %get.Length Get method for the dependent property 'Length'
%             
%             rd  = obj.RData;
%             val = rd(end);
%         end
    end
    
    methods % construction
        function obj = Stick(varargin)
            %Stick Constructor for the class 'Stick'
            
            %Pass it on
            obj@awi.model.Component(varargin{:});    
            
            %Extend property groups
            obj.addPropertyGroup('Appearance', ...
                'LineColor', 'Colour applied to stick model'    , ...
                'LineStyle', 'Line style applied to stick model', ...
                'LineWidth', 'Width of the line');
            
        end
    end
    
    methods % visualisation
        
        function hg = drawElement(obj, ht, tag)
            
            %Caller supplied tag explicitly?
            if nargin < 3
                
                %No - apply a tag indicating that we are just drawing some lines
                tag = 'Sticks';
                
            end
            
            %Start with base class
            hg = drawElement@awi.model.Component(obj, ht);
            
            %Get x,y,z data carefully (e.g. so user can edit content bit by bit)
            [xd, yd, zd] = xyzdata(obj);
            
            %Get additional arguments
            args = { ...
                'Color'          , obj.LineColor       ; ...
                'LineStyle'      , obj.LineStyle       ; ...
                'LineWidth'      , obj.LineWidth       ; ...
                'Marker'         , obj.MarkerStyle     ; ...
                'MarkerSize'     , obj.MarkerSize      ; ...
                'MarkerFaceColor', obj.MarkerFaceColor ; ...
                'MarkerEdgeColor', obj.MarkerEdgeColor ; ...
                'Tag'            , tag}';
            
            %Eliminate any empties
            args(:, cellfun(@isempty, args(2,:))) = [];
            
            %Draw the line
            hg{end+1} = line('Parent', ht, ...
                'XData', xd, ...
                'YData', yd, ...
                'ZData', zd, ...
                args{:});
            
        end
        
    end
    
    methods (Static) % getLineLength
        
        function r = getLineLength(x, y, z, r0)
            %getLineLength Calculates the straight line distance between
            %the points (x,y,z) and returns a vector of distances starting
            %at 'r0'.
            
            if nargin < 4
                r0 = 0;
            end
            
            %Simply the cumulative sum of the hypotenuse lengths!
            dXYZ = diff(x).^2 + diff(y).^2 + diff(z).^2;
            if iscolumn(dXYZ)
                r = cumsum([r0 ; sqrt(dXYZ)]);
            else
                r = cumsum([r0, sqrt(dXYZ)]);
            end
            
        end
        
    end
    
    methods % s2pos, xyzdata        
        function pos = s2pos(obj, s, flag)
            %s2pos Interpolates the coordinates of the stick using the
            %options defined by:
            %   - 's'    : Normalised position along the parent span
            %   - 'flag' : Defines the span vector along which the
            %              interpolation will take place.
            
            %Grab x,y,z data
            [xd, yd, zd] = xyzdata(obj);
            
            %Check to see if the line data for the parent object has been
            %defined - If not then return (0, 0, 0) and move on
            if isempty(xd)
                pos = [0, 0, 0];
                return
            end
            
            switch flag
                case {'X', 'Y', 'Z'}
                    %Define eta - Depends on offset axis definition
                    switch flag
                        case 'X'
                            eta = abs(xd ./ max(abs(xd)));
                        case 'Y'
                            eta = abs(yd ./ max(abs(yd)));
                        case 'Z'
                            eta = abs(zd ./ max(abs(zd)));
                    end
                    %Interpolate!
                    pos(:, 1) = interp1(eta, xd, s);
                    pos(:, 2) = interp1(eta, yd, s);
                    pos(:, 3) = interp1(eta, zd, s);
                case 'R'
                    %Invoke get method once!
                    rd = obj.RData;
                    %Define eta
                    eta = rd ./ max(rd);
                    %Interpolate!
                    r  = interp1(eta, rd, s);
                    pos(:, 1) = interp1(rd, xd, r);
                    pos(:, 2) = interp1(rd, yd, r);
                    pos(:, 3) = interp1(rd, zd, r);
            end
            
        end
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

