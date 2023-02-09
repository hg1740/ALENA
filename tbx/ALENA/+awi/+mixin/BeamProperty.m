classdef BeamProperty < matlab.mixin.SetGet
    %BeamProperty Describes a quantity that can be defined over a beam.
    %
    % A beam property has the following characterisitics:
    %
    %   1. 'Quantity' - This is the name of the beam property. It must be a
    %   row vector of characters.
    %
    %   2. 'Value' - This is the value of the beam property at some
    %   position along the beam. It must be a numeric row vector.
    %
    %   3. 'Distribution' - This is the points along the beam where 'Value'
    %   is defined.
    %
    %   4. 'Variation' - This describes the variation of the value of the
    %   beam property between the defined sample points. It must be a valid
    %   method for the function 'interp1'.
    %
    %   5. 'Type'  - 
    %   
    %   6. 'Name' - 
    %
    %   7. 'Description' - 
    %
    %   8. 'Units'
    %
    % Method list:
    %
    %   * 'getEta' - This method obtains the unique points along the beam
    %   by combining all 'Distributions' from an array of
    %   awi.mixin.BeamProperty objects. This method can be called on a
    %   scalar or object array, if called on a scalar object then it will
    %   return the value of 'Distribution'.
    %
    %   * 'getBPV' - This method interpolates a beam property at user-
    %   specified sample points using user-specified arguments for the
    %   function 'interp1'.
    % 
    % TODO
    %   - Add a 'DistributionType' property for handling continuous or
    %   discontinous data sets. e.g. Cover thicknesses, etc.
    
    %Labelling/Naming the beam property
    properties %Name, Description, Units
        %'Name' is the 'pretty' string presented to the user in the view
        Name = '';
        %'Description' is a character array describing the quantity
        Description  = '';
        %'Units' are the the units that the quantity is defined in
        Units = '';
        %'Label' is the text that will appear in the y-axis label
        Label = '';
    end
    
    %Beam property parameters
    properties %Quantity, Value, Distribution, Variation, AxisFlag, Type 
        %'Quantity' is the name of the beam property
        Quantity     = '';
        %'Value' is the value of the beam property at the points defined in
        %'Distribution'.
        Value        = [nan , nan];
        %'Distribution' is the sample points along the beam. It is
        %non-dimensionsal between 0 and 1.
        Distribution = [0 , 1];
        %'Variation' is how the value of the beam property varies between
        %the sample points.
        Variation    = 'linear';
        %'AxisFlag' is the axis along which the value of 'eta' is defined.
        %It be one of the 3 global axes (X,Y,Z) or the local axis 'R'.
        AxisFlag = 'R';
        %'Type' denotes a sub-group of beam properties
        Type = '';
        %Interpolated value - Placeholder required by 'Beamable'.
        Interpolated 
    end
    
    %Beam geometry
    properties %x, y, z
        %Beam x-coordinates
        x
        %Beam y-coordinates
        y
        %Beam z-coordinates
        z
    end
    
    properties (Dependent) %etaX, etaY, etaZ, etaR, BeamGeom, SegLength, r, Length
        %Normalised coordinates along the global X-axis
        etaX
        %Normalised coordinates along the global Y-axis
        etaY
        %Normalised coordinates along the global Z-axis
        etaZ
        %Normalised coordinates along the beam axis
        etaR
        %Combination of the (x,y,z) data of the beam. Order of
        %concatenation is [x ; y ; z].
        BeamGeom
        %Length of each straight line section of the beam
        SegLength
        %Straight line distance along the beam through all the (x,y,z)
        %points defined in 'X'.
        r
        %Total length of the beam when measuring along the beam axis
        Length
    end
    
    methods % set / get 
        function set.Quantity(obj, val)      %set.Quantity     
            %set.Quantity Set method for the property 'Quantity'.
            %
            % 'Quantity' is the name of the beam property.
            validateattributes(val, {'char'}, {'row'}, class(obj), 'Quantity');
            obj.Quantity = val;
        end
        function set.Value(obj, val)         %set.Value        
            %set.Value Set method for the property 'Value'.
            %
            % 'Value' is the value of the beam property at the points along
            % the beam defined at 'Distribution'.            
            
            %TODO - Remove this once Phil has fixed the bug with
            %'validateattributes' being caught by a try-catch statement
            if isnumeric(val) && iscolumn(val)
               val = val'; 
            end
            
            if isscalar(val)
                val = repmat(val, size(obj.Distribution)); %#ok<MCSUP>
            end
            
            %Verify it is a matrix
            validateattributes(val, {'numeric'}, {'2d'}, class(obj), 'Value');   
            
            %Verify dimensions
            [d1, d2] = size(val);
            assert(d1 >= 1, ['The number of rows in ''Value'' must ', ...
                'be greater than or equal to 1']);
            assert(d2 >= 2, ['The number of columns in ''Value'' '  , ...
                'must be greater than or equal to 2']);
            
            %Assign the value!
            obj.Value = val;
        end
        function set.Distribution(obj, val)  %set.Distribution 
            %set.Distribution Set method for the property 'Distribution'.
            %
            % 'Distribution' is the sample points along the beam which must
            % be non-dmensionalised between 0 and 1.
            
            val = obj.validateDistribution(obj, val, 'Distribution');
            if isscalar(val)
                error(['The distribution of the beam property ''%s'' ', ...
                    'must have at least 2 points along the beam ', ...
                    'otherwise it is not possible to interpolate the ', ...
                    'beam properties.'], obj.Quantity); %#ok<MCSUP>
            end
            obj.Distribution = val;
        end
        function set.Variation(obj, val)     %set.Variation    
            %set.Variation Set method for the property 'Variation'.
            %
            % 'Variation' is the how the value of the beam property varies
            %between the sample points.
            
            obj.validateVariation(obj, val);
            obj.Variation = val;
        end
        function set.AxisFlag(obj, val)      %set.AxisFlag     
           %set.AxisFlag Set method for the property 'AxisFlag'.
           %
           % 'AxisFlag' is the axis along which the beam property is
           % defined. It can be one of the 3 global axes or the local 'R'
           % axis.
           %validatestring(val, {'X', 'Y', 'Z', 'R'}, class(obj), 'AxisFlag');
           obj.validateAxisFlag(obj, val);
           obj.AxisFlag = upper(val);
        end
        function set.Type(obj, val)          %set.Type         
            %set.Type Set method for the property 'Type'.
            %
            % 'Type' is the name of a subgroup of beam properties.
            validateattributes(val, {'char'}, {'row'}, class(obj), 'Type');
            obj.Type = val;
        end    
        function set.x(obj, val)             %set.x            
           %set.x Set method for the property 'x'.
           %
           % 'x' must be a numeric row vector of real values.
           validateattributes(val, {'numeric'}, {'row', 'real',  ...
               'finite', 'nonnan', 'nonempty'}, class(obj), 'x');
           obj.x = val;
        end
        function set.y(obj, val)             %set.y            
           %set.y Set method for the property 'x'.
           %
           % 'y' must be a numeric row vector of real values.
           validateattributes(val, {'numeric'}, {'row', 'real',  ...
               'finite', 'nonnan', 'nonempty'}, class(obj), 'y');
           obj.y = val;
        end
        function set.z(obj, val)             %set.z            
           %set.z Set method for the property 'z'.
           %
           % 'z' must be a numeric row vector of real values.
           validateattributes(val, {'numeric'}, {'row', 'real',  ...
               'finite', 'nonnan', 'nonempty'}, class(obj), 'z');
           obj.z = val;
        end
        function val = get.Name(obj)         %get.Name         
            val = obj.Name;
            if isempty(val)
                val = obj.Quantity;
            end
        end
        function val = get.Label(obj)        %get.Label        
            val = obj.Label;
            if isempty(val)
                val = [obj.Name, ' [', obj.Units, ']'];
            end
        end
        function val = get.Units(obj)        %get.Units        
            val = obj.Units;
            if isempty(val)
                val = '-';
            end
        end
        function val = get.Type(obj)         %get.Type         
            %get.Type Get method for the property 'Type'.
            %   
            %   Returns the string "- Unassigned Beam Property -" if the
            %   value has not already been set by the user during the
            %   initiation of the beam property.
            
            val = obj.Type;
            if isempty(val)
                val = '- Unassigned Beam Property -';
            end
        end
        function val = get.BeamGeom(obj)     %get.BeamGeom     
            %get.BeamGeom Get method for the dependent property 'BeamGeom'.
            
            %Retrive data in the appropriate form
            [xd, yd, zd] = xyzdata(obj);            
            %Concatente
            val = [xd ; yd ; zd];
        end
        function val = get.SegLength(obj)    %get.SegLength    
            %get.SegLength Get method for the dependent property
            %'SegLength'.
            
            %Each section is a straight line so the length of each section
            %is just the sum of the squared differences.
            val = sqrt(sum(diff(obj.BeamGeom, [], 2).^2));
        end
        function val = get.r(obj)            %get.r            
            %get.R Get method for the interpolated property 'R'.
            % 
            % 'R' is the straight line distance between all (x,y,z) points
            % defined in 'obj.X'.
            
            %Simply the cumulative sum of the hypotenuse lengths!
            sl  = obj.SegLength;
            val = cumsum([0, sl]);
            
        end
        function val = get.Length(obj)       %get.Length       
           %get.Length Get method for the dependent property 'Length'.
           val = obj.r(end);
        end
        function val = get.etaX(obj)         %get.etaX         
            if isempty(obj.x)
                val = [];
                return
            end
            val = obj.x ./ obj.x(end);
        end
        function val = get.etaY(obj)         %get.etaY         
            if isempty(obj.y)
                val = [];
                return
            end
            val = obj.y ./ obj.y(end);
        end
        function val = get.etaZ(obj)         %get.etaZ         
            if isempty(obj.z)
                val = [];
                return
            end
            val = obj.z ./ obj.z(end);
        end
        function val = get.etaR(obj)         %get.etaR         
           if isempty(obj.BeamGeom)
               val = [];
               return
           end
           val = obj.r ./ obj.Length;
        end
    end
    
    methods % constructor
        function obj = BeamProperty(varargin)
            %BeamProperty Constructor for the 'awi.mixin.BeamProperty'
            %class.            

            %Grab any properties from 'varargin' and assign to object
            %   - Assume varargin in in parameter name/value format
            prps  = properties(obj);
            names = varargin(1 : 2 : end);
            value = varargin(2 : 2 : end);
            idx   = ismember(names, prps);
            set(obj, names(idx), value(idx));
            
        end
    end
        
    methods (Sealed) % getEta, getBPV, selectBeamProp
        
        % sampling beam properties
        function val = getEta(obj, type, varargin)
            %getEta Shortcut method for obtaining the finest possible 
            %distrbution of points along the notional beam by combining all
            %'Distribution' values for each beam property in the object
            %array 'obj' which are of type 'type'.
            %
            % If the beam property type is not passed in as the second
            % argument then Eta is unique combinaton of all 'eta' points
            % along the beam.
            
            if iscolumn(obj)
                obj = obj';
            end
            if nargin < 2
               type = 'all'; 
            end
            
            %Parse inputs
            p = inputParser;
            addRequired(p , 'obj'     , @(x)validateattributes(x, {'awi.mixin.BeamProperty'}, {'row'}));
            addRequired(p , 'type'    , @(x)i_validatetype(x, obj)); 
            addParameter(p, 'AxisFlag', 'R'  , @(x)any(validatestring(x, {'X', 'Y', 'Z', 'R'})));
            parse(p, obj, type, varargin{:});
            
            %Extra validation for the type input
            function i_validatetype(type, obj)
                
                types = unique({obj.Type});
                names = {obj.Quantity};
                
                %Create anonymous function for parsing 'type'
                valid_fcn = @(x)(validatestring(x,  [{'all'}, types, names], class(obj), 'type'));
                
                if iscell(type) 
                    idx_ = cellfun(@(x) valid_fcn(x), type, 'Unif', false); %#ok<NASGU>
                else
                    valid_fcn(type);
                end                
                
            end
            
            %Grab inputs from parser
            type   = p.Results.type;
            axflag = p.Results.AxisFlag;
            
            %Down-select beam properties
            if any(strcmpi(type, 'all'))
                %Use 'Distribution' of ALL beam properties
                BP = obj;
            else
                %Down-select the beam properties based on beam property type
                BP  = selectBeamProp(obj, type);
            end
            
            %Do we need to adjust any distributions?
            idx = ismember({BP.AxisFlag}, axflag);
            
            distr       = cell(size(BP));
            distr(idx)  = {BP(idx).Distribution};            
            distr(~idx) = arrayfun(@(bp) switchEtaAxis(bp, axflag), BP(~idx), 'Unif', false);
            
            %All unique values - Use a tolerance to avoid floating point
            %error
            tol = 1e-10;
            %val = unique([horzcat(distr{:}), [0, 1]]);
            val = uniquetol([horzcat(distr{:}), [0, 1]], tol);
            
        end
        function val = getBPV(obj, Quantity, samplePoints, method, varargin)
            %getBPV Retrieves the Beam Property Value (BPV) for a given 
            %'Quantitity'. The BPV is the value of the 'Quantity' field 
            %which is interpolated at 'samplePoints' using the 
            %interpolation method 'method'.
            %
            % Inputs 
            %   - 'obj'
            %     Array of 'awi.mixin.BeamProperty' objects.
            %
            %   - 'samplePoints'
            %     Row vector of sample points in the range [0, 1].
            %
            %   - 'method'  
            %     Valid interpolation method for use with 'interp1'.
            %
            %   - 'varargin' 
            %     Variable argument list of valid parameters for 'interp1'.
                        
            %Collect the valid BeamProperty names
            validBeamProp = {obj.Quantity};
            
            %Check that 'Quantity' is valid
            validatestring(Quantity, validBeamProp, class(obj), ...
                'the requested beam property');
            
            %Grab BeamProperty object related for the requested 'Quantity'
            idx = ismember(validBeamProp, Quantity);
            BP  = obj(idx);
            
            %Has the user supplied a value of 'AxisFlag'?
            idx = ismember(varargin, 'AxisFlag');
            if ~any(idx)
                axFlag = BP.AxisFlag;
            else
                ind = find(idx == true);
                axFlag = varargin{ind + 1};
                %Strip from 'varargin' as we no longer need it
                varargin([ind, ind+1]) = []; 
            end
            
            %Have the sample points been passed to the function?
            if nargin < 3
                %If not, default to the finest possible distribution of
                %points along the beam for this beam property type
                samplePoints = getEta(obj, BP.Type, 'AxisFlag', axFlag);
            end 
            
            %Has an interpolation method been provided?
            if nargin < 4
               method =  BP.Variation;
            end
            if isempty(method)
                method = BP.Variation;
            end
            
            %Have any interpolation parameters been defined?
            if isempty(varargin)
                %Extrapolate data as standard
                varargin{1} = 'extrap';
            end

            %Check 'Value' & 'Distribution' have the same no. of elements
            if numel(BP.Distribution) ~= size(BP.Value, 2)
                if strcmp(BP.Name, 'RMatrix')
                    val = nan(size(BP.Distribution));
                    return
                end
                error('The number of points in the distribution and the value must be the same');
            end
            
            %Check that the 'Distribution' does not contain any duplicates
            %   - TODO : Need to verify if this works for matrix of
            %   'Distribution' values.
            [xv, ia, ~] = unique(BP.Distribution);
            if numel(xv) ~= numel(BP.Distribution)
                %Grab the data for the unique points
                yv = BP.Value(ia);                
                %Get index number of duplicate values
                index = (1 : numel(BP.Distribution))'; %Index of all points
                idx_d = ~ismember(index, ia);
                ind_d = index(idx_d);
                %For each duplicate point -> Use the average of the
                %duplicate values
                repVal = zeros(1, numel(ind_d));
                for iD = 1 : numel(ind_d)
                    %Find all occurences of the duplicate points
                    idx = ismember(BP.Distribution, BP.Distribution(ind_d(iD)));
                    %ind = find(BP.Distribution == BP.Distribution(ind_d(iD)));
                    %Compute the average and store it
                    repVal(iD) = mean(BP.Value(idx));
                end
                %Index into the new distribition to determine the mapping
                %for the new 'smeared' data
                idx = ismember(xv, BP.Distribution(ind_d));
                %Assign the 'smeared' data to the 'yv' data set
                yv(idx) = repVal;
            else
                %No duplicate points so just pass the values through to the
                %interpolation!
                yv = BP.Value;
            end           
            
            %Perform the interpolation
            %   - Tranpose to allow vectorised interpolation of multiple
            %     data sets.
            val = interp1(xv', yv', samplePoints', method, varargin{:})';
        end
        
        % down-selecting beam properties
        function BP = selectBeamProp(obj, type)
            %selectBeamProp Returns all the beam properties of Type 'type'.
            
            %Gather all beam types
            allTypes  = {obj.Type};
                        
            %Index all beam properties belonging to this 'type'.
            idx = ismember(allTypes, type);
                        
            %Check to see if a beam property was actually selected
            if ~any(idx)
                %Maybe the user requested a beam property by name?
                allNames = {obj.Quantity};
                idx      = ismember(allNames, type);
            end
            
            %Select the BeamProperty object
            BP  = obj(idx);
            
        end
        
    end
 
    methods (Access = private) %xyzdata, switchEtaAxis 
        function [xd, yd, zd] = xyzdata(obj)
            
            %In principle it's as simple as this
            xd = obj.x;
            yd = obj.y;
            zd = obj.z;
            
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
        function newEta = switchEtaAxis(obj, newEtaAxis)
            %switchEtaAxis Returns the values of 'eta' along the axis
            %specified by 'newEtaAxis'.
            %
            % The coordinates of the beam are defined at points (xi,yi,zi)
            % and the length along the beam is given by ri.
            %
            % The beam property is defined at j points ('distr') along an 
            % axis specified by 'obj.AxisFlag'.
            %
            % The new eta vector is calculated by finding the (x,y,z,r)
            % coordinates at all the j points ('distr') along the beam and
            % then normalising the coordinates by the appropriate
            % i-coordinate as defined by 'newEtaAxis'.
            %
            % TODO - Update this so that is uses vector instead of trig
            
            %Check for the case where the beam has no geometry
            if obj.Length == 0 
                newEta = obj.Distribution;
                return
            end
            
            %Grab (x,y,z) & r data
            [xi, yi, zi] = xyzdata(obj);
            ri = obj.r;
                        
            %Segment lengths?
            dx_i = diff(xi, [], 2);
            dy_i = diff(yi, [], 2);
            dz_i = diff(zi, [], 2);
            dr_i = sqrt(dx_i.^2 + dy_i.^2 + dz_i.^2);
            
            %Orientation angles for each segment
            alpha = atan2(dz_i, dx_i); %#ok<NASGU>
            beta  = atan2(dy_i, dx_i);
            gamma = atan2(dz_i, dy_i);
            
            %What is the current distribution?
            distr = obj.Distribution .* obj.(lower(obj.AxisFlag))(end);
            
            %Calculate the (x,y,z,r) coordinates at all j points
            switch obj.AxisFlag
                case 'X'
                    L  = [0, cumsum(dx_i)]; %Cumulative length
                    dl = dx_i;              %Increment
                    %Find the segment index and the fraction of each point
                    %along that segment
                    [f, n] = i_findSegmentFraction(distr, L, dl);
                    %Convert the fraction to a length along the beam
                    dx_j = f .* dl(n);
                    %Calculate (x,y,z,r)
                    dr_j = dx_j ./ cos(beta(n));
                    r_j  = dr_j + ri(n);
                    x_j  = dx_j + xi(n);
                    y_j  = sin(beta(n))  .* r_j + yi(n);
                    z_j  = sin(gamma(n)) .* r_j + zi(n);
                case 'Y'
                    L  = [0, cumsum(dy_i)]; %Cumulative length
                    dl = dy_i;              %Increment
                    %Find the segment index and the fraction of each point
                    %along that segment
                    [f, n] = i_findSegmentFraction(distr, L, dl);
                    %Convert the fraction to a length along the beam
                    dy_j = f .* dl(n);
                    %Calculate (x,y,z,r)
                    dr_j = dy_j ./ sin(beta(n));
                    r_j  = dr_j + ri(n);
                    x_j  = cos(beta(n)) .* dr_j + xi(n);
                    y_j  = dy_j + yi(n);
                    z_j  = cos(gamma(n)) .* dr_j + zi(n);
                case 'Z'
                    L  = [0, cumsum(dz_i)]; %Cumulative length
                    dl = dz_i;              %Increment
                    %Find the segment index and the fraction of each point
                    %along that segment
                    [f, n] = i_findSegmentFraction(distr, L, dl);
                    %Convert the fraction to a length along the beam
                    dz_j = f .* dl(n);
                    %Calculate (x,y,z,r)
                    dr_j = dz_j ./ sin(gamma(n));
                    r_j  = dr_j + ri(n);
                    x_j  = cos(beta(n)) .* dr_j + xi(n);
                    y_j  = cos(beta(n)) .* dr_j + yi(n);
                    z_j  = dz_j + zi(n);
                case 'R'
                    L  = [0, cumsum(dr_i)]; %Cumulative length
                    dl = dr_i;              %Increment
                    %Find the segment index and the fraction of each point
                    %along that segment
                    [f, n] =  i_findSegmentFraction(distr, L, dl);
                    %Convert the fraction to a length along the beam
                    dr_j = f .* dl(n);
                    %Calculate (x,y,z)
                    x_j = cos(beta(n))  .* dr_j + xi(n);
                    y_j = sin(beta(n))  .* dr_j + yi(n);
                    z_j = sin(gamma(n)) .* dr_j + zi(n);
                    r_j = dr_j + ri(n);
            end
            
            %Normalise the data along the appropriate axis to obtain the
            %new eta distribution
            switch newEtaAxis
                case 'X'
                    newEta = x_j ./ xi(end);
                case 'Y'
                    newEta = y_j ./ yi(end);
                case 'Z'
                    newEta = z_j ./ zi(end);
                case 'R'
                    newEta = r_j ./ ri(end);
            end
            
            %Remove nan terms
            newEta(isnan(newEta)) = [];
            
            %TODO - Debug code to understand why we are getting an eta
            %value > 1
            newEta(newEta > 1) = 1;
            
            function [f, n] =  i_findSegmentFraction(distr, L, dl)
                %i_findSegmentFraction Finds the position of a point along
                %a segment as a fraction of the segment length.
                %
                % TODO - Vectorise this!!
                
                n = zeros(size(distr));
                f = zeros(size(distr));
                for i = 1 : numel(distr)
                    if distr(i) == L(1)
                        f(i) = 0;
                        n(i) = 1;
                    elseif distr(i) == L(end)
                        f(i) = 1;
                        n(i) = numel(L) - 1;
                    else
                        n(i) = find(distr(i) > L, 1, 'last');
                        f(i) = 1 - (L(n(i) + 1) - distr(i)) ./ dl(n(i));
                    end
                end                
            end
        end
    end
    
    methods (Static) %validating properties 
        function val = validateDistribution(obj, val, tok)
            %validateDistribution Checks that the distribution is formatted
            %correctly. 
            %
            % i.e. Numeric vector with all values between 0 <= val <= 1.
            
            if isnumeric(val) && iscolumn(val)
               val = val'; 
            end            
            validateattributes(val, {'numeric'}, {'row', 'nonnegative', ...
                '<=', 1, 'nonempty', 'finite', 'real', 'nonnan'}, class(obj), tok);            
        end
        function val = validateVariation(obj, val)
           %validateVariation Checks the value of 'Variation' is correct.
           %
           % 'Variation' must be one of the valid variation options
           % provided by 'interp1'. 
           %
           % See also, INTERP1
           
            validMethod = {'linear', 'nearest', 'next', 'previous', ...
                'spline', 'pchip', 'cubic', 'v5cubic'};
            val = validatestring(val, validMethod, class(obj), 'Variation');
            
        end
        function validateAxisFlag(obj, val)
            %validateAxisFlag Checks the value of 'AxisFlag' is correct.
            %
            % This is a static method so that it can be accessed by the
            % 'awi.mixin.Beamable' method that listens to changes in the
            % BeamPropertyObject values.
            %
            % 'AxisFlag' is the axis along which the beam property is
            % defined. It can be one of the 3 global axes or the local 'R'
            % axis.
            validatestring(val, {'X', 'Y', 'Z', 'R'}, class(obj), 'AxisFlag');
        end
    end
    
end

