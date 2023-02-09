classdef MassDistribution < mvc.view.LayeredDrawing
    %MassDistribution Shows a detailed hierarchical break-down of the mass
    %distribution of a given model.
    %
    % Mass has three components: 
    %   - Point masses  : These are described by the 'awi.model.PointMass'
    %                     class.
    %   - Element Mass  : Dependent on the type of element. For beam (1D)
    %                     elements this is found by integrating the 
    %                     Non-Structural Mass (NSM).
    %   - Material Mass : Found by integrating the product of the material
    %                     density distribution and the element 
    %                     cross-sectional area.


    properties (Access = protected) %Graphics handles
        %Handle to the graphics objects representing the masses
        hMass
        %Handle to the uicontrol object for changing the scale factor
        hScaleFactor
        %Handle to the uicontrol object for toggling the scaling of masses
        hScaleMass
    end
    
    properties %ScaleFactor
        %Factor for scaling the radius of the spheres
        ScaleFactor = 1;
        %Logical flag for toggling the scaling of masses
        ScaleMass   = true;
    end
    
    methods % set / get
        function set.ScaleFactor(obj, val) %set.ScaleFactor
            %set.ScaleFactor Set method for the property 'ScaleFactor'.
            %
            %   - 'ScaleFactor' must be a scalar numeric.
            %
            
            %Escape route
            if isempty(val)
                return
            end
            
            %Validate
            validateattributes(val, {'numeric'}, {'scalar', 'nonnan', ...
                'finite', 'real', 'nonempty'}, class(obj), 'ScaleFactor');
            obj.ScaleFactor = val;
            
        end
        function set.ScaleMass(obj, val)   %set.ScaleMass
            %set.ScaleMass Set method for the property 'ScaleMass'.
            %
            %   - 'ScaleMass' must be a scalar logical.
            %
            
            %Escape route
            if isempty(val)
                return
            end
            
            %Validate
            validateattributes(val, {'logical'}, {'scalar'}, ...
                class(obj), 'ScaleMass');
            obj.ScaleMass = val;
            
        end
    end
    
    methods % construction        
        function obj = MassDistribution(varargin)
            
            %Call superclass constructor
            obj@mvc.view.LayeredDrawing(varargin{:});
           
            %Prevent the user from changing which layer of the model is 
            %visible or invisible
            obj.hDrawLayers.Enable  = 'off';
            obj.hDrawLayers.Visible = 'off';
            
            %Do NOT link this view to the GUI selection as it is expensive
            %to generate
            obj.SelectionLinked        = false;
            obj.hSelectionLinked.Value = false;
            
            %Add an option for the user to change the scale-factor
            uicontrol('Parent', obj.hGrid, 'Style', 'text', ...
                'String'             , [10, ' Scale Factor'], ...
                'HorizontalAlignment', 'left');            
            obj.hScaleFactor = uicontrol('Parent', obj.hGrid, ...
                'Style', 'edit', 'String', '1', ...
                'Callback', @obj.cbUpdateMassRadius);  
            
            %Add an option for toggling scaling the mass values
            obj.hScaleMass = uicontrol('Parent', obj.hGrid, ...
                'Style'   , 'check'     , ...
                'String'  , 'Scale Mass', ...
                'Value'   , true        , ...
                'Callback', @obj.cbUpdateMassRadius);
            
            %Set appearance            
            obj.hGrid.Widths = -1;
            obj.hGrid.Heights(end - 2 : end) = 25;
            
        end
    end
    
    methods (Access = protected) % updating the view
        function [L, ha] = update(obj)
            %update Updates the view.
            
            %Force sticks to be visible 
            if ~any(ismember(obj.DrawLayers, 'Sticks'))
                obj.DrawLayers{end + 1} = 'Sticks';
            end
            
            %Delete any existing mass graphical objects
            if ~isempty(obj.hMass)
                delete(obj.hMass);
            end
            
            %Start with base-class
            [L, ha] = update@mvc.view.LayeredDrawing(obj);
            
            %Find content that is currently visible
            hc = findobj(ha, 'Visible', 'on');
            
            %Get tagged content
            tag = get(hc, {'Tag'}); 
            
            %Only interested if tag non-empty
            b = cellfun(@isempty, tag);
            tag(b) = [];
            hc(b)  = [];
            
            %Identify tags that do NOT match what we're looking for
            b = ~ismember(tag, 'Sticks');
            
            %Hide that content
            set(hc(b), 'Visible', 'off');      
            
            %Turn off the 'Selected' property
            setProperty(L, 'Selected', 'off', ha);
            
            %Grab the 'awi.model.PointMass' objects & get the CoG location
            switch obj.DrawOption
                
                case 'All'
                    
                    %Easy - Just get to the top of the tree and flatlist
                    frmwrk = obj.Selection.Root;
                    lst    = flatlist(frmwrk);
                    bMass  = arrayfun(@(x) isa(x, 'awi.model.PointMass'), lst);
                    Masses = lst(bMass);
                    
                    %Calculate the CoG (including all children)
                    cg = frmwrk.Aircraft.CoG;
                    
                    %Construct the axes title
                    titstr = sprintf('CoG location for the model is (%#.4g, %#.4g, %#.4g)', ...
                        cg(1), cg(2), cg(3));
                    
                case 'Selection and children'
                    
                    %Easy - Just flatlist everything from here down
                    lst    = flatlist(obj.Selection);
                    bMass  = arrayfun(@(x) isa(x, 'awi.model.PointMass'), lst);
                    Masses = lst(bMass);
                    
                    %Calculate the CoG (including all children)
                    cg = obj.Selection.CoG;
                    
                    %Construct the axes title
                    titstr = sprintf(['CoG location for the %s '    , ...
                        '(including all child mass elements) is '   , ...
                        '(%#.4g, %#.4g, %#.4g)'], obj.Selection.Name, ...
                        cg(1), cg(2), cg(3));
                    
                case 'Selection only'
                    
                    %Look for any (rogue) 'awi.model.PointMass' objects
                    %   - Should really be collected by
                    %   'awi.model.PointMasses' objects
                    bMass  = arrayfun(@(x) isa(x, 'awi.model.PointMass'), ...
                        obj.Selection.Children);                    
                    Masses = obj.Selection.Children(bMass);
                    
                    %Check for any collected 'awi.model.PointMass' objects
                    bMasses = arrayfun(@(x) isa(x, 'awi.model.PointMasses'), ...
                        obj.Selection.Children);
                    
                    %Combine
                    Masses = [Masses ; vertcat(obj.Selection.Children(bMasses).Children)];
                    
                    %Calculate the CoG (including all children)
                    cg = obj.Selection.ThisCoG;
                    
                    %Construct the axes title
                    titstr = sprintf(['CoG location for the %s '    , ...
                        '(excluding all child mass elements) is '   , ...
                        '(%.4g, %.4g, %.4g)'], obj.Selection.Name, ...
                        cg(1), cg(2), cg(3));
                    
                otherwise
                    error('Bad selection');
            end
            
            %Update axes title
            ha.Title.String = titstr;
            
            %Escape route
            if isempty(Masses)
                return
            end
            
            %How many masses?
            nMasses = numel(Masses);
            
            %Get the position of the point masses
            pos = vertcat(Masses.AbsPosition);
            
            %Get the mass of the point masses
            M = [Masses.Mass];
            
            %Get the coordinates of a unit sphere
            [X, Y, Z] = sphere;
            
            %Repeat the cooordinates by 'nMasses'
            X = repmat(X, [1, 1, nMasses]);
            Y = repmat(Y, [1, 1, nMasses]);
            Z = repmat(Z, [1, 1, nMasses]);
            
            %Normalise the masses?
            if obj.ScaleMass
                %Normalise the mass values and combine with the scaling factor
                %to link the appearance of the spheres to their mass
                R = (M ./ max(M)) .* obj.ScaleFactor;
            else
                R = repmat(obj.ScaleFactor, [1, nMasses]);
            end            
            
            %Ensure correct size for matrix-wise operation!
            R = permute(R, [3, 1, 2]);  
            
            %Scale the sphere(s)
            X = X .* R;
            Y = Y .* R;
            Z = Z .* R;
            
            %Translate the sphere(s)
            X = X + permute(pos(:, 1), [3, 2, 1]);
            Y = Y + permute(pos(:, 2), [3, 2, 1]);
            Z = Z + permute(pos(:, 3), [3, 2, 1]);
            
            %Plot the masses as spheres
            set(ha, 'NextPlot', 'add'); %Make sure we don't replace the data
            hM = arrayfun(@(i) surf(ha, X(:, :, i), Y(:, :, i), Z(:, :, i), ...
                'FaceColor', 'b', 'FaceAlpha', 0.75), 1 : nMasses, 'Unif', false);  
            
            %Store the handles
            obj.hMass = horzcat(hM{:});
            
            %Stash the mass of each sphere in the user data
            set(obj.hMass, {'UserData'}, num2cell(M)');
                       
        end        
        function cbUpdateMassRadius(obj, ~, ~)
            %cbUpdateMassRadius Updates the appearance of the spheres. 
            %
            % If the user has selected to scale the mass radii then the 
            % radius of the spheres will be scaled by the value in the 
            % 'editbox' 'obj.hScaleFactor'.
            
            %Update the scale factor with the value from the uicontrol
            obj.ScaleFactor = str2double(obj.hScaleFactor.String);
            obj.ScaleMass   = logical(obj.hScaleMass.Value);
            
            %Grab the coordinates and user data from the spheres
            data = get(obj.hMass, {'XData', 'YData', 'ZData', 'UserData'});
            
            %Return original data to matrix format for indexing
            X = cat(3, data{:, 1});
            Y = cat(3, data{:, 2});
            Z = cat(3, data{:, 3});
            M = cat(2, data{:, 4});
            
            %Grab the dimensions
            [d1, d2, d3] = size(X);
            
            %Find origin of each mass
            x_ = permute(arrayfun(@(i) mean(mean(X(:, :, i))), 1 : d3), [1, 3, 2]);
            y_ = permute(arrayfun(@(i) mean(mean(Y(:, :, i))), 1 : d3), [1, 3, 2]);
            z_ = permute(arrayfun(@(i) mean(mean(Z(:, :, i))), 1 : d3), [1, 3, 2]);
            
            %Remove the translation before scaling
            X = X - x_;
            Y = Y - y_;
            Z = Z - z_;
            
            %Find the radius of each sphere
            xR = permute(arrayfun(@(i) max(max(X(:, :, i))), 1 : d3), [1, 3, 2]);
            yR = permute(arrayfun(@(i) max(max(Y(:, :, i))), 1 : d3), [1, 3, 2]);
            zR = permute(arrayfun(@(i) max(max(Z(:, :, i))), 1 : d3), [1, 3, 2]);
            
            %Return the spheres to a unit sphere
            X = X ./ xR;
            Y = Y ./ yR;
            Z = Z ./ zR;
            
            %Normalise the masses?
            if obj.ScaleMass
                %Normalise the mass values and combine with the scaling factor
                %to link the appearance of the spheres to their mass
                R = (M ./ max(M)) .* obj.ScaleFactor;
            else
                R = repmat(obj.ScaleFactor, [1, d3]);
            end            
                        
            %Ensure correct size for matrix-wise operation!
            R = permute(R, [3, 1, 2]);  
            
            %Scale the sphere(s)
            X = X .* R;
            Y = Y .* R;
            Z = Z .* R;
            
            %Return back to their original positions
            X = X + x_;
            Y = Y + y_;
            Z = Z + z_;
            
            %Convert 3D matrix into cell arrays
            X = squeeze(mat2cell(X, d1, d2, ones(1, d3)));
            Y = squeeze(mat2cell(Y, d1, d2, ones(1, d3)));
            Z = squeeze(mat2cell(Z, d1, d2, ones(1, d3)));
            
            %Apply back to the graphics objects
            set(obj.hMass, {'XData'}, X);
            set(obj.hMass, {'YData'}, Y);
            set(obj.hMass, {'ZData'}, Z);
            
        end
    end
    
end

