classdef BoxBeam < awi.model.Component & awi.mixin.BeamPropertyObject
    %BoxBeam Defines a cross-section which provides beam geometric
    %properties. 
    %
    % * This class in analogous to the MSC.Nastran PBARL bulk data type.
    %
    % * The following box-beam types are available:
    %   - Box: Describes a  
    %   - 
    %
    % TODO - If I subclass this from 'awi.model.CrossSection' then I will
    % have difficulties adding dynamic properties to the parent
    % 'awi.model.Beam' object.
    
    %Box-Beam type
    properties
        %Box Beam idealisation
        BoxType = 'SymmetricBox';
    end
    
    %Parameters for calculating component stresses
    properties
        %Method for calculating the shear stress
        ShearStressMethod = 'thin-wall';
        %Transverse Load Factor (Symmetric Box only)
        TransverseLoadFactor = 0.5;
    end
    
    %Properties related to the various box-beam types
    properties (Constant = true)
        %Allowable box types
        ValidBoxTypes   = {'SymmetricBox', 'Box'};
        %Allowable shear stress methods 
        ValidShrStrMethod = {'thin-wall', 'simplified'};
        %Properties relating to the various box properties
        BoxProperties   = { ...
            {'Width', 'Height', 'CoverThickness', 'SparThickness'}, ...
            {'Width', 'Height', 'TopCoverThickness', 'BottomCoverThicknesss', 'FrontSparThickness', 'RearSparThickness'}};
        %Shorthand labels for the box properties
        BoxPropLabel = { ...
            {'W', 'H', 'TC' , 'TS'}, ...
            {'W', 'H', 'TTC', 'TBC', 'TFS', 'TRS'}};
        %Methods for calculating the geometric properties
        BoxMethods       = {@calcSymmetricBoxProps, @calcBoxProps};
        %Methods for calculating the component stress
        StressMethods    = {@calcSymmetricBoxStress, @calcBoxStress};
        DirStressMethods = {};
        ShrStressMethods = {};
        %Equivalent names of the cross-sections in MSC.Nastran
        NastranBoxCode   = {'BOX', 'BOX1'};
        %Equivelent names of the cross-section properties in MSC.Nastran
        NastranPropCode  = { ...
            {'DIM1', 'DIM2', 'DIM3', 'DIM4'}, ...
            {'DIM1', 'DIM2', 'DIM3', 'DIM4', 'DIM6', 'DIM5'}}
        %Equivalent names of the MSC.Nastran beam properties that can be
        %calculated
        BeamPropMap = { ...
            {'A' , 'I1' , 'I2' , 'J'  , 'NSM', 'NSI' ; ...  %Nastran Beam 
            'Abb', 'Izz', 'Ixx', 'Jbb', 'NSM', 'NSI'}, ...  %AWI Box-Beam 
            {'A' , 'I1' , 'I2' , 'I12', 'J'  , 'NSM', 'NSI' ; ...   %Nastran Beam
            'Abb', 'Izz', 'Ixx', 'Ixz', 'Jbb', 'NSM', 'NSI'}};      %AWI Box-Beam
    end
    
    %Beam geometric properties
    properties
        %Area of the cross-section which is effective in carrying axial
        %stress
        Abb = nan;
        %Second moment of area in the 'xx' plane
        Ixx = nan;
        %Second moment of area in the 'zz' plane
        Izz = nan;
        %Cross second moment of area in the 'xz' plane
        Ixz = nan;
        %Polar moment of area which is effective in torsion
        Jbb = nan;
        %Non-structural mass
        NSM = nan;
        %Non-structural inertia
        NSI = nan;
        %Location of the neutral axis w.r.t. the BoxBeam origin
        xNA = nan;
        zNA = nan;
        %Location of the shear centre w.r.t the BoxBeam origin
        xSC = nan;
        zSC = nan;
        %Location of the centre of mass w.r.t the BoxBeam origin
        xCM = nan;
        zCM = nan;
    end
    
    %Box-Beam geometry
    properties
       %Height of the section in the local z-direction [m]
       Height = nan;
       %Width of the section in the local x-direction [m]
       Width = nan;
       %Thickness of the top and bottom covers [m]
       CoverThickness = nan;
       %Thickess of the front and rear spar [m]
       SparThickness = nan;
    end
   
    %Generic Box-Beam geometry
    properties
        %(x,z) coordinate of the local BoxBeam origin w.r.t. the beam axis
        SectionOrigin = [0 ; 0];
        %(x,z) coordinates of the corner points
        %   - Coordinates are defined in the clockwise direction starting
        %     at the bottom left corner of the box model. (i.e. the leading
        %     edge/front spar).
        CornerCoords  = [nan ; nan];
        %Thickness of each piecewise linear component
        Thickness     = nan;
    end
    
    %Generic Box-Beam geometry
    properties 
        %(x,z) coordinates of the corner points (in a continuous loop).
        SectionCoords
        %Length of each piecewise-linear component
        EdgeLength
        %Slenderness ratio of each piecewise linear component
        BuponT
        %x-coordinates of the centroid of each piecewise linear component
        xC
        %x-coordinates of the centroid of each piecewise linear component
        zC
    end
    
    %Box Properties
    properties (Dependent)
        CurrentPropNames
        CurrentPropLabel
        CurrentNasBoxCode
        CurrentNasPropCode        
        CurrentNasPropNames
    end
    
    methods % set / get
        function set.BoxType(obj, val)              %set.BoxType
            %set.BoxType Set method for the property 'BoxType'.
            %
            % 'BoxType' must be one of the prescribed valid box types as
            % defined by 'obj.ValidBoxType'.
            
            validatestring(val, obj.ValidBoxTypes, class(obj), 'BoxType');
            obj.BoxType = val;
        end
        function set.ShearStressMethod(obj, val)    %set.ShearStressMethod
            %set.ShearStressMethod Set method for the property
            %'ShearStressMethod'.
            %
            % 'ShearStressMethod' must be one of the methods defined by
            % 'obj.ValidShrStrMethod'.
            
            validatestring(val, obj.ValidShrStrMethod, class(obj), 'ShearStressMethod');
            obj.ShearStressMethod = val;
            
        end
        function set.TransverseLoadFactor(obj, val) %set.TransverseLoadFactor
            %set.TransverseLoadFactor(obj, val) Set method for the property
            %'TransverseLoadFactor'. 
            %
            % 'TransverseLoadFactor' must be a scalar numeric between [0:1]
            
            validateattributes(val, {'numeric'}, {'scalar', '<=', ...
                1, '>=', 0}, class(obj), 'TransverseLoadFactor');
            obj.TransverseLoadFactor = val;
        end
        function set.SectionOrigin(obj, val)        %set.SectionOrigin
            %set.SectionOrigin Set method for the property 'SectionOrigin'.
            %
            % 'SectionOrigin' must be a 2x1 numeric vector
            
            validateattributes(val, {'numeric'}, {'size', [2, 1], ...
                'real', 'finite', 'nonnan'}, class(obj), 'SectionOrigin');
            obj.SectionOrigin = val;
        end
        function set.CornerCoords(obj, val)         %set.CornerCoords
            %set.CornerCoords Set method for the property 'CornerCoords'.
            
            validateattributes(val, {'numeric'}, {'nrows', 2, ...
                'nonnan', 'finite', 'real'}, class(obj), 'CornerCoords');
            obj.CornerCoords = val;
        end
        function set.Thickness(obj, val)            %set.Thickness
            %set.Thickness Set method for the property 'Thickness'.
            
            validateattributes(val, {'numeric'}, {'nonnegative', ...
                'row', 'finite', 'real', 'nonnan'}, class(obj), 'Thickness');
            obj.Thickness = val;
        end
        function val = get.SectionCoords(obj)       %get.SectionCoords
            %get.SectionCoords Get method for the property 'SectionCoords'.
            
            val = obj.SectionCoords;
            
            %'SectionCoords' is the same as 'CornerCoords' but the final
            %coordinates are repeated to create a continuous loop.
            if isempty(val)
                val = [obj.CornerCoords, obj.CornerCoords(:, 1)];
            end
            
        end
        function val = get.CurrentPropNames(obj)    %get.CurrentPropNames
            idx = ismember(obj.ValidBoxTypes, obj.BoxType);
            val = obj.BoxProperties{idx};
        end
        function val = get.CurrentPropLabel(obj)    %get.CurrentPropLabel
            idx = ismember(obj.ValidBoxTypes, obj.BoxType);
            val = obj.BoxPropLabel{idx};
        end
        function val = get.CurrentNasBoxCode(obj)   %get.CurrentNasBoxCode
            idx = ismember(obj.ValidBoxTypes, obj.BoxType);
            val = obj.NastranBoxCode{idx};
        end
        function val = get.CurrentNasPropCode(obj)  %get.CurrentNasPropCode
            idx = ismember(obj.ValidBoxTypes, obj.BoxType);
            val = obj.NastranPropCode{idx};
        end
        function val = get.CurrentNasPropNames(obj) %get.CurrentNasPropNames
            idx = ismember(obj.ValidBoxTypes, obj.BoxType);
            val = obj.BeamPropMap{idx}(1, :); 
        end
    end
    
    methods % constructor
        function obj = BoxBeam(varargin)
            %BoxBeam Constructor for the 'awi.model.BoxBeam' class.
            %
            % Actions performed:
            %   - Call superclass constructor
            
            %Pass it on
            obj@awi.model.Component(varargin{:});
            
        end
    end
    
    methods % visualisation
        function hg = drawElement(obj, ha, mode)
            %drawElement Object-specific visualisation method.
            
            if nargin < 2 || isempty(ha)
                hF  = figure('Name', 'Box-Beam Cross-Section');
                ha  = axes('Parent', hF, 'NextPlot', 'add');
            end
            if nargin < 3
                mode = 'local';
            end
            validatestring(mode, {'local', 'global'}, class(obj(1)), 'mode');
            
            switch mode 
                case 'global'
                    
                    %Grab coordinates
                    [x, z] = getLocalCoords(obj);         
                    
                    %Generate 3D data
                    y = repmat([obj.BeamEta], [size(x, 1), 1]) .*  obj(1).BeamHandle.Span;      
                    
                    %Plot box edges
                    hg = plot3(ha, x, y, z);
                    
                case 'local'
                    
                    %Grab coordinates
                    [x, z] = getSectionCoords(obj);
                    [x, z] = addLocalOffsets(obj, x, z);
                    
                    %Plot box edges
                    hg = plot(ha, x', z');
                    
                    %Format axes
                    xlabel(ha, 'X [m]');
                    zlabel(ha, 'Z [m]');
                    axis(ha, 'equal');
                    
                    %Plot the centroids
                    xC_ = vertcat(obj.xC)';
                    zC_ = vertcat(obj.zC)';
                    [xC_, zC_] = addLocalOffsets(obj, xC_, zC_);
                    plot(ha, xC_, zC_, ...
                        'LineStyle'      , 'None', ...
                        'Marker'         , 'o'   , ...
                        'MarkerFaceColor', 'r'   , ...
                        'MarkerEdgeColor', 'k')
                    
            end
           
            function [x, z] = getLocalCoords(obj)
                %getLocalCoords Retrieves the (x, z) coordinates in the
                %local box coordinate system.
                
                %Grab data
                [x, z] = getCoords(obj);                
                x = [x, x(:, 1)]';
                z = [z, z(:, 1)]';
                
                %Add local offsets
                origin = [obj.SectionOrigin];
                x = x + origin(1, :);
                z = z + origin(2, :);
                
            end
            
        end
    end
    
    methods % public facing methods for running analysis
        %Splitting corner coords into (x, z) components
        function [x, z] = getCoords(obj)
            %getCoords Splits the 'CornerCoords' property into seperate
            %(x,z) vectors.
            
            x = [];
            z = [];
            
            %Grab data - remove empties
            coords = {obj.CornerCoords};
            coords(cellfun(@isempty, coords)) = [];
            
            if isempty(coords)
                return
            end
            
            %Split into (x,z) components
            x = cellfun(@(x) x(1, :), coords, 'Unif', false);
            z = cellfun(@(x) x(2, :), coords, 'Unif', false);
            x = vertcat(x{:});
            z = vertcat(z{:});
            
        end
        function [x, z] = getSectionCoords(obj)
            %getCoords Splits the 'SectionCoords' property into seperate
            %(x,z) vectors.
            
            [x, z] = getCoords(obj);
            if isempty(x) || isempty(z)
                return
            end
            x = [x, x(:, 1)];
            z = [z, z(:, 1)];
            
        end
        %Calculating geometric properties
        function getGeometricProps(obj)
            %getGeometricProps Calculates the values of the geometric
            %propertes for the chosen box-beam type.
            
            %Check the user has requested something sensible
            validatestring(obj(1).BoxType, obj(1).ValidBoxTypes);
            
            %Check all the required properties have been defined
            checkBoxBeamProperties(obj);
            
            %Get handle to the method which will update the properties for
            %the chosen box type
            fn = obj(1).BoxMethods{ismember(obj(1).ValidBoxTypes, obj(1).BoxType)};
            
            %Do it
            fn(obj);
            
            %Calculate any 'private' properties
            updateSectionProperties(obj);

        end
        %Calculating component stress (direct & shear)
        function [sig, tau, vm] = getComponentStress(obj, Mxx, Mzz, AxF, Fx, Fz, Trq)
            %getComponentStress Calculates the stress in each component of
            %the box beam model.
            
            nObj = numel(obj);            
            flag = 'both';
            
            if nargin < 2   %Check for correct number of loads
                warning('No load data provided to the ''awi.model.BoxBeam'' objects. Unable to calculate component stress.');
                return
            elseif nargin == 2  
                %Allow user to input a matrix of loads
                loads = Mxx;
                sz    = size(loads);
                idx   = sz == nObj;
                %Check format of input data
                if numel(sz) ~=2    %Correct dimensions?
                    i_printWarning;
                    return
                end
                if ~any(idx)        %Correct size?
                    i_printWarning;
                    return
                elseif find(idx, 1, 'first') == 1
                    loads = loads';
                end
                nLoad = sz(~idx);
                if nLoad == 3       %Correct no. loads?
                    Mxx  = loads(1, :);
                    Mzz  = loads(2, :);
                    AxF  = loads(3, :);
                    flag = 'direct';
                elseif nLoad == 6
                    Mxx  = loads(1, :);
                    Mzz  = loads(2, :);
                    AxF  = loads(3, :);
                    Fx   = loads(4, :);
                    Fz   = loads(5, :);
                    Trq  = loads(6, :);
                    flag = 'both';
                else
                    i_printWarning;
                    return
                end
                clear loads
            elseif nargin < 5   
                %Direct stress only 
                flag = 'direct';
            end
                   
            %Sensible defaults
            sig = [];
            tau = [];
            vm  = [];
            
            %Check the user has requested something sensible
            validatestring(obj(1).BoxType, obj(1).ValidBoxTypes);
            
            %Check all the required properties have been defined
            checkBoxBeamProperties(obj);
            
            %Get handle to the method which will calculate the stress for
            %the chosen box type
            switch flag
                case 'direct'
                    fn = obj(1).DirStressMethods{ismember(obj(1).ValidBoxTypes, obj(1).BoxType)};                    
                    fn(obj,  Mxx, Mzz, AxF, Fx, Fz, Trq);    
                case 'both'
                    fn = obj(1).StressMethods{ismember(obj(1).ValidBoxTypes, obj(1).BoxType)};                    
                    [sig, tau, vm] = fn(obj, Mxx, Mzz, AxF, Fx, Fz, Trq);                    
            end
            
            function i_printWarning
                warning('The correct number of loads must be provided. Unable to calculate component stresses.');
            end
            
        end
    end
    
    methods (Access = private) % internal methods which do the actual analysis
        function [x, z] = addLocalOffsets(obj, x, z)
            %addLocalOffsets Add the local section offsets to each (x, z)
            %coordinate.            
            
            nObj = numel(obj);            
            sz   = size(x);
                        
            %Add local offsets
            origin = [obj.SectionOrigin];
            dx     = origin(1, :);
            dz     = origin(2, :);
            if find(sz == nObj, 1, 'first') == 1
                dx = dx';
                dz = dz';
            end            
            x = x + dx;
            z = z + dz;
                
        end
        %Calculating geometric properties
        function calcSymmetricBoxProps(obj)
            %calcSymmetricBoxProps Calculates the geometric properties
            %for the 'SymmetricBox' box-beam type.
            %
            % The 'SymmetricBox' has is symmetric about the (x,z) axes.
            % Therefore, the cross-second moment of area (Ixz) is zero.
            % i.e. There is no coupling between in-plane and out of plane
            % coupling.
            %
            % The equations for the A, Ixx, Izz, Ixz & J values are taken
            % from the MSC.Nastran Reference Manual PBEAML entry for the
            % bulk data card (TYPE = BOX).

            %Validate input
            if iscolumn(obj)
                obj = obj';
            end
            validateattributes(obj, {'awi.model.BoxBeam'}, {'row'}, ...
                class(obj), 'obj');
            
            %Can we get the density value?
            if isempty(obj(1).BeamHandle)
                %Assume the density value
                %   - TODO: Update this so the user can pass in the density
                %   value of we should just error out at this point.
                rho = 2810;
            else
                if isempty(obj(1).BeamHandle.Material)
                    rho = nan;
                else
                    rho = obj(1).BeamHandle.Material(1).Rho;
                end
            end
            
            %Shorthand variables
            zrs = zeros(size(obj));
            h   = [obj.Height];
            w   = [obj.Width];
            tc  = [obj.CoverThickness];
            ts  = [obj.SparThickness];
            hi  = h - (2 .* tc); %height of the inner square
            wi  = w - (2 .* ts); %width of the inner square
            
            %Calculate the properties
            a   = (h .* w) - (hi .* wi);
            izz = ((h .* w.^3) ./ 12) - ((hi .* wi.^3) ./ 12);
            ixx = ((w .* h.^3) ./ 12) - ((wi .* hi.^3) ./ 12);
            ixz = zeros(size(h));
            j   = ((2 .* tc .* ts) .* (w - ts).^2 .* (h - tc).^2) ./ ...
                (w .* ts + h .* tc - tc.^2 - ts.^2);
            nsm = rho .* (w .* h - wi .* hi);
            
            %Break down NSI calculation as it is a bit more involved...
            %   - Cover
            mc = rho .* tc .* wi;
            rc = (h - tc) ./ 2;
            Ic = (mc ./ 12) .* (wi.^2 + tc.^2) + mc .* rc.^2;
            %   - Spar
            ms = rho .* ts .* h;
            rs = (w - ts) ./ 2;
            Is = (ms ./ 12) .* (h.^2 + ts.^2) + ms .* rs.^2;
            nsi = 2 .* (Ic + Is);
            
%             nsi = 2 .* rho .* ( ...
%                 (ts .* h .* (((h.^2 + ts.^2)/12) + ((w - ts).^2/4))) + ...  %spar
%                 (tc .* wi .* (((wi.^2 + tc.^2)/12) + ((h - tc)/4))));       %cover
            
            %Assign the data to the objects
            set(obj, {'Abb'}, num2cell(a)');
            set(obj, {'Ixx'}, num2cell(ixx)');
            set(obj, {'Izz'}, num2cell(izz)');
            set(obj, {'Ixz'}, num2cell(ixz)');
            set(obj, {'Jbb'}, num2cell(j)');
            set(obj, {'NSM'}, num2cell(nsm)');
            set(obj, {'NSI'}, num2cell(nsi)');
            
            %As the box is symmetric the offsets are in the middle
            xOff = num2cell(w ./ 2)';
            zOff = num2cell(h ./ 2)';
            set(obj, {'xSC'}, xOff);
            set(obj, {'zSC'}, zOff);
            set(obj, {'xNA'}, xOff);
            set(obj, {'zNA'}, zOff);
            set(obj, {'xCM'}, xOff);
            set(obj, {'zCM'}, zOff);
            
            %Assign the thickness values
            tVec = [ts ; tc ; ts ; tc]';
            set(obj, {'Thickness'}, num2cell(tVec, 2));
            
            %Calculate the corner coordinates
            x = [zrs ; zrs ; w ; w  ]';
            z = [zrs ; h   ; h ; zrs]';
            coords = arrayfun(@(i) [x(i, :) ; z(i, :)], 1 : numel(obj), 'Unif', false);
            set(obj, {'CornerCoords'}, coords');
            
        end
        function calcBoxProps(obj)
            %calcBoxProps Calculates the geometric properties for the
            %'Box' box-beam type.
            
            error('Update code');
        end        
        function checkBoxBeamProperties(obj)
            %checkBoxBeamProperties Checks that the necessary box-beam
            %properties have been defined.
            
            assert(numel(unique({obj.BoxType})) == 1, ['The value of ', ...
                '''BoxType'' must be the same for all BoxBeam objects.']);
            
            %Get the properties
            boxprops = obj(1).BoxProperties{ismember(obj(1).ValidBoxTypes, obj(1).BoxType)};
            val      = get(obj, boxprops);
            val      = cell2mat(val); %TODO - This will fall over if any of the box properties are non-numeric
            
            %Check for nan
            if numel(obj) == 1
                idx_nan = isnan(val);
            else
                idx_nan = any(isnan(val));
            end
            
            %Assert that there are no nan entries
            assert(~any(idx_nan), ['Error calculating the box-beam ', ...
                'properties for BoxType = ''%s''\n\nThe following ' , ...
                'box-beam properties are undefined:\n\n\t%s\n']     , ...
                obj(1).BoxType, strjoin(boxprops(idx_nan), ', '));
            
        end
        %Calculating component stress (Symmetric Box)
        function [sig, tau, vm] = calcSymmetricBoxStress(obj, Mxx, Mzz, AxF, Fx, Fz, Trq)
            %calcSymmetricBoxStress Calculates the shear and direct stress
            %in each component of the BoxBeam section.

            sig = calcSymmetricBoxDirectStress(obj, Mxx, Mzz, AxF);
            tau = calcSymmetricBoxShearStress(obj , Fx , Fz , Trq);
            
            vm = sqrt(sig .^ 2 + 3 .* tau .^ 2);
            
        end
        function sig = calcSymmetricBoxDirectStress(obj, Mxx, Mzz, AxF)
            %calcSymmetricBoxDirectStress Calculates the direct stress in
            %each BoxBeam component for the 'SymmetricBox' idealisation.
            
            %Grab section properties
            A   = [obj.Abb]; %#ok<*PROPLC>
            Ixx = [obj.Ixx];
            Izz = [obj.Izz];
            Ixz = [obj.Ixz];
            MoA = Ixx .* Izz - Ixz .^ 2;
            x   = [obj.Width]  ./ 2;
            z   = [obj.Height] ./ 2;
            
            %Calculate direct stress
            sig = ...
                (((Mzz .* Ixx - Mxx .* Ixz) ./ MoA) .* x) + ...
                (((Mxx .* Izz - Mzz .* Ixz) ./ MoA) .* z) + (AxF ./ A);
            
        end
        function tau = calcSymmetricBoxShearStress(obj, Fx, Fz, Trq)
            %calcSymmetricBoxShearStress Calculates the shear stress in
            %each BoxBeam component for the 'SymmetricBox' idealisation.
            
            %Grab section properties
            tC = [obj.CoverThickness];
            tS = [obj.SparThickness];
            h  = [obj.Height];
            w  = [obj.Width];
            k  = [obj.TransverseLoadFactor];
                   
            %Area enclosed by the thin walls
            Ai = (w - tS) .* (h - tC);
            
            tau = zeros(2, numel(obj));
            fTau = @(t, A, h, w, Trq, F1, F2, k) ...
                (abs(Trq) ./ (2 .* A .* t))   + ...
                (abs(F1 ) ./ (2 .* h .* t))   + ...
                (k .* abs(F2) ./ (2 .* w .* t));
            
            %Shear stress in cover & spar
            tau(1, :) = fTau(tC, Ai, h, w, Trq, Fx, Fz, k);
            tau(2, :) = fTau(tS, Ai, h, w, Trq, Fz, Fx, k);
            
        end
    end
    
    methods (Access = private) % generic box-beam methods
        function updateSectionProperties(obj)
            %updateSectionProperties Updates the following properties of
            %the BoxBeam section.
            %   - Segment length (b)
            %   - Segment slenderness ratio (b/t)
            %   - Segment centroid (xB, zB)
            %   - First Moment of Area (Qxx, Qzz)
            %   - Neutral Axis location (xNA, zNA)
            %   - Second Moment of Area (Ixx, Izz, Ixz)
            %   - Polar Moment of Area (J)
            
            getLineLength(obj);     %   - b
            getBuponT(obj);         %   - b/t
            getEdgeCentroid(obj);   %   - (xB, zB)
            
            
        end
        function getLineLength(obj)
            %getLineLength Calculates the length of each piecewise linear
            %component.
            
            [x, z] = getSectionCoords(obj);            
            r = sqrt(diff(x, [], 2) .^ 2 + diff(z, [], 2) .^2);            
            set(obj, {'EdgeLength'}, num2cell(r, 2));
            
        end
        function getBuponT(obj)
            %getBuponT Calculates the b/t ratio for each edge.
            
            b = {obj.EdgeLength};
            b = vertcat(b{:});
            t = {obj.Thickness};
            t = vertcat(t{:});    
            
            bUponT = b ./ t;            
            
            set(obj, {'BuponT'}, num2cell(bUponT, 2));
            
        end
        function getEdgeCentroid(obj)
            %getEdgeCentroid Calculates the centroid coordinates of each
            %piecewise lienar segment.
            
            [x, z] = getSectionCoords(obj);
            
            %xB_N = x_N + (x_N+1 - x_N) ./ 2 
            xB = x(:, 1 : end - 1) + (diff(x, [], 2) ./ 2); 
            zB = z(:, 1 : end - 1) + (diff(z, [], 2) ./ 2); 
            
            set(obj, {'xC'}, num2cell(xB, 2));
            set(obj, {'zC'}, num2cell(zB, 2));
            
        end
    end
    
end

