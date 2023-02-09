classdef LoadDistribution < matlab.mixin.SetGet
    %LoadDistribution Describes a distribution of loads (forces & moments)
    %and enforced motion (displacements, velocities & accelerations)
    %along the length of a beam.
    %
    % Detailed Description:
    %
    %  * This object describes a SINGLE distribution of applied loads and
    %  enforced motions along a beam (i.e. one set of nondimensional
    %  points). If the user desires a different distribution of points then
    %  they can modify the 'EtaDistribution' property or create a new
    %  object.
    %
    %   * The number of points in 'EtaDistribution' (nP) must match the
    %   number of columns in the applied loads and enforced motions.
    %   
    %   * Time domain loading can be specified by adding addition pages to
    %   the applied loads and enforced motions - e.g. expand the 2D
    %   matrices into 3D matrices, however,the number of pages must equal
    %   the number of time domain points (nT)
    %
    %   * If you are unsure if the dimensions of your load/motion
    %   distribution are correct then call the method ''parse'' on the
    %   object to check the dimensions.
    %       i.e. >> L = awi.loads.LoadDistribution >> parse(L)
    %
    % TODO - Add offsets for loads!!
    
    %Definition of loads/motions along the beam and across time
    properties
        %Nondimensionsal position of all loads and enforced motions along
        %the beam. -> numel(EtaDistribution) = nP
        EtaDistribution
        %Denotes the beam axis which 'EtaDistribution' is defined along
        AxisFlag = 'R';
        %Time vector for all loads and enforced motions
        % -> numel(TimeVector) = nT
        TimeVector = 0;
    end
    
    %Applied loads
    properties
        %Loads applied at a point along the beam [6, nP, nT]
        PointLoads
        %Coordinate system in which the point loads are defined
        PointLoadCoordSys = 'local';
        %Behaviour of the load
        PointLoadBehaviour = 'non-follower';
        %Loads defined at a point along the beam but assumed to vary
        %linearly from one point to the next [6, nP, nT]
        DistributedLoads
        %Coordinate system in which the distributed loads are defined
        DistrLoadCoordSys = 'local';
    end
    
    %Enforced motion
    properties
        %Enforced displacement along the beam [6, nP, nT]
        Displacement
        %Enforced velocity along the beam [6, nP, nT]
        Velocity
        %Enforced acceleration along the beam [6, nP, nT]
        Acceleration
    end
    
    %Visualisation
    properties
        %Maximum length of the load vector
        MaximumVectorLength = 1;
    end
    
    %Handling relationship with 'awi.model.Beam' objects
    properties (SetAccess = private)
        %Reference to the beam handle which this object belongs to.
        BeamHandle
    end
    
    %Internal properties
    properties (Dependent, Hidden)
        %Coordinates of the applied load w.r.t the root of the beam
        DistributionCoords
    end
    
    %Helper properties
    properties (Constant)
        %Valid coordinate system tokens
        ValidCoordSys = {'local', 'global'};
    end
    
    methods % set / get
        function set.EtaDistribution(obj, val)     %set.EtaDistribution
            %set.EtaDistribution Set method for the property
            %'EtaDistribution'.
            
            if isempty(val)%User trying to reset?
                obj.EtaDistribution = [];
                return
            end
            validateattributes(val, {'numeric'}, {'row', 'nonnegative', ...
                '<=', 1, 'nonnan', 'real', 'increasing'}, class(obj), 'EtaDistribution');
            obj.EtaDistribution = val;
        end
        function set.TimeVector(obj, val)          %set.TimeVector
            %set.TimeVector Set method for the property 'TimeVector'.
            
            if isempty(val) %User trying to reset?
                val = 0;
            end
            validateattributes(val, {'numeric'}, {'row', 'nonnegative', ...
                'nonnan', 'real'}, class(obj), 'EtaDistribution');
            obj.TimeVector = val;
        end
        function set.AxisFlag(obj, val)            %set.AxisFlag
            %set.AxisFlag Set method for the property 'AxisFlag'.
            
            awi.mixin.BeamProperty.validateAxisFlag(obj, val);
            obj.AxisFlag = upper(val);
        end
        function set.PointLoads(obj, val)          %set.PointLoads
            %set.PointLoads Set method for the property 'PointLoads'.
            
            val = validateLoadDistribution(obj, val, 'PointLoads');
            obj.PointLoads = val;
        end
        function set.PointLoadCoordSys(obj, val)   %set.PointLoadCoordSys
            %set.PointLoadCoordSys Set method for the property
            %'PointLoadCoordSys'.
            
            validatestring(val, obj.ValidCoordSys, class(obj), 'PointLoadCoordSys');
            obj.PointLoadCoordSys = lower(val);
        end
        function set.PointLoadBehaviour(obj, val)  %set.PointLoadBehaviour
            %set.PointLoadBehaviour Set method for the property
            %'PointLoadBehaviour'.
            val = lower(val);
            validatestring(val, {'follower', 'non-follower'}, ...
                class(obj), 'PointLoadBehaviour');
            obj.PointLoadBehaviour = val;
        end
        function set.DistributedLoads(obj, val)    %set.DistributedLoads
            %set.DistributedLoads Set method for the property 'DistributedLoads'.
            
            val = validateLoadDistribution(obj, val, 'DistributedLoads');
            obj.DistributedLoads = val;
        end
        function set.DistrLoadCoordSys(obj, val)   %set.PointLoadCoordSys
            %set.DistrLoadCoordSys Set method for the property
            %'DistrLoadCoordSys'.
            
            validatestring(val, obj.ValidCoordSys, class(obj), 'DistrLoadCoordSys');
            obj.DistrLoadCoordSys = lower(val);
        end
        function set.Displacement(obj, val)        %set.Displacement
            %set.Displacement Set method for the property 'Displacement'.
            
            val = validateLoadDistribution(obj, val, 'Displacement');
            obj.Displacement = val;
        end
        function set.Velocity(obj, val)            %set.Velocity
            %set.Velocity Set method for the property 'Velocity'.
            
            val = validateLoadDistribution(obj, val, 'Velocity');
            obj.Velocity = val;
        end
        function set.Acceleration(obj, val)        %set.Acceleration
            %set.Acceleration Set method for the property 'Acceleration'.
            
            val = validateLoadDistribution(obj, val, 'Acceleration');
            obj.Acceleration = val;
        end
        function set.BeamHandle(obj, val)          %set.BeamHandle
            %set.BeamHandle Set method for the property 'BeamHandle'.
            
            %User may be detaching the load distribution from the beam
            if isempty(val)
                val = [];
            else
                validateattributes(val, {'awi.model.Beam'}, {'scalar'}, ...
                    class(obj), 'BeamHandle');
            end
            obj.BeamHandle = val;
        end
        function set.MaximumVectorLength(obj, val) %set.MaximumVectorLength
            %set.MaximumVectorLength Set method for the property
            %'MaximumVectorLength'.
            
            validateattributes(val, {'numeric'}, {'scalar', 'real', ...
                'positive'}, class(obj), 'MaximumVectorLength');
            obj.MaximumVectorLength = val;
        end
        function val = get.DistributionCoords(obj) %get.DistributionCoords
            %DistributionCoords Get method for the dependent property
            %'DistributionCoords'.
            
            %Pass it on
            val = getDistributionCoords(obj);
        end
    end
    
    methods % relating the load distribution to a beam
        function assignLoadToBeam(obj, BeamObj)
            %assignLoadToBeam Assigns the load distribution to a beam
            %object and vice-versa.
            
            set(obj, 'BeamHandle', BeamObj);
            if isempty(BeamObj) %Escape route
                %User is detaching load
                return
            end
            assignLoadToBeam(BeamObj, obj);
        end        
    end
     
    methods % visualisation
        function hg = drawElement(obj, ht, tag)
            %drawElement Object-specific method for visualising the
            %load/motion distribution along the beam.
            %
            % Cannot proceed if the load distribution is not assigned to a
            % beam.
            
            hg = [];
            
            if nargin < 3
                tag = 'Load/Motion Distribution';
            end
            
            if numel(obj) > 1
                error('Update code to allow vectorised plotting of object arrays');
            end
            
            %Check that the dimensions of the properties are appropriate
            parse(obj);
            
            if isempty(obj.BeamHandle) %Escape route
                return
            end
            
            nLoads = numel(obj.EtaDistribution);            
            
            %Get the (x,y,z) coordinates of the distribution
            xyz = getDistributionCoords(obj);
            
            %Point loads
            pf = obj.PointLoads(1 : 3, :, 1);   %forces - TODO Update for all time histories
            pm = obj.PointLoads(4 : 6, :, 1);   %moments
            %   - Normalise
            pf = obj.MaximumVectorLength .* pf ./ max(max(abs(pf)));
            %   - Transform into local frame?
            if strcmpi(obj.PointLoadCoordSys, 'local')
                %Get orientation at load points
                rotMat = interp1(obj.BeamHandle.Orientation_eta, ...
                    obj.BeamHandle.RMatrix_', obj.EtaDistribution, 'previous');
                %Transform into local frame
                for i = 1 : nLoads
                   pf(:, i) = reshape(rotMat(i, :), [3, 3]) * pf(:, i);                    
                end
            end
            %   - Place the loads at the correct point along the beam
            pf = xyz + pf;
            %   - Generate line coordinates
            [pfX, pfY, pfZ] = i_padWithNaN(xyz, pf);
            %   - Plot it
            hg{1} = line('Parent', ht, ...
                'XData', pfX, ...
                'YData', pfY, ...
                'ZData', pfZ, ...
                'Marker', 'none', ...
                'LineStyle', '-', ...
                'Color', 'm', ...
                'Tag', 'Point Loads');
            
            %Distributed loads...
            %   - TODO Plot with a patch!!
            
            function [x, y, z] = i_padWithNaN(xyzLoad0, xyzLoad)
                %i_padWithNaN Pads the coordinates with nan terms to allow
                %vectorised plotting.
                
                nan_ = nan(1, size(xyzLoad0, 2));
                
                %Pad with nan terms
                x = [xyzLoad0(1, :) ; xyzLoad(1, :) ; nan_];
                y = [xyzLoad0(2, :) ; xyzLoad(2, :) ; nan_];
                z = [xyzLoad0(3, :) ; xyzLoad(3, :) ; nan_];
                
                %Return column vector
                x = x(:);
                y = y(:);
                z = z(:);
                
            end
                        
        end
    end
    
    methods (Access = private) % helper functions
        function xyz = getDistributionCoords(obj)
            %getDistributionCoords Returns the (x,y,z) coordinates of the
            %load/motion distribution w.r.t the root of the assocaited beam
            %object.
            
            %Sensible default
            xyz = [];
            
            %Cannot proceed if there is not beam associated with the object
            if isempty(obj.BeamHandle)
                return
            end
            
            %Grab the beam geometry data
            [xd, yd, zd] = xyzdata(obj.BeamHandle);            
            if isempty(xd) || isempty(yd) || isempty(zd) %Escape route
                return
            end
            
            %Interpolate the data based on which global axis the eta
            %distribution is defined along
            switch obj.AxisFlag
                case 'X'    %Global X
                    eta = xd ./ xd(end);
                case 'Y'    %Global Y
                    eta = yd ./ yd(end);
                case 'Z'    %Global Z
                    eta = zd ./ zd(end);
                case 'R'    %Straight-line distance along beam                   
                    r   = awi.model.Stick.getLineLength(xd, yd, zd);
                    eta = r ./ r(end);
            end
            xyz = interp1(eta, [xd ; yd ; zd]', obj.EtaDistribution)';
            
        end
    end
    
    methods % validation
        function parse(obj)
            %parse Checks that the load/motion definition matches the
            %spatial and temporal definitions prescribed by
            %'EtaDistribution' and 'TimeVector'.
            %
            % The following properties are checked:
            %   - PointLoads
            %   - DistributedLoads
            %   - Displacement
            %   - Velocity
            %   - Acceleration
            %
            % If one of the properties is empty then it is assumed that the
            % user is not utlising that property and it is skipped.
            
            if numel(obj) > 1
                arrayfun(@parse, obj)
                return
            end
            
            if isempty(obj.EtaDistribution)
                warning(['The spatial distribution (''EtaDistribution'') '    , ...
                    'has not been defined. Specify the spanwise distribution ', ...
                    'of the applied loads and enforced motions then proceed.']);
                return
            end
            
            %Required dimensions
            nP = numel(obj.EtaDistribution);
            nT = numel(obj.TimeVector);
            
            %Properties to be parsed
            prp = {'PointLoads', 'DistributedLoads', 'Displacement', 'Velocity', 'Acceleration'};
            
            %Required dimensions
            tok = {'Load/Motion', 'EtaDistribution', 'TimeVector'};
            dim = [6, nP, nT];
            str = arrayfun(@(i) sprintf(' (dim%i / numel = %i)', i, dim(i)), ...
                1 : numel(dim), 'Unif', false);
            
            %Check sizes
            for iP = 1 : numel(prp)
                if isempty(obj.(prp{iP}))
                    continue
                end
                [d1, d2, d3]  = size(obj.(prp{iP}));
                idx = ([d1, d2, d3] == dim);
                msg = strjoin(strcat(tok(~idx), str(~idx)), ', ');
                assert(all(idx), sprintf(['Expected the property ''%s'' ', ...
                    'to be of size [6, %i, %i]. Update the following ', ...
                    'dimensions before proceeding:\n\t- %s\n\n'], ...
                    prp{iP}, nP, nT, msg));
            end
            
        end
    end   
    
    methods (Access = private) % validation
        function val = validateLoadDistribution(obj, val, tok)
            %validateLoadDistribution Checks that the attributes of 'val',
            %which describes the property 'tok', match the required format.
            %The required attributtes of a generic load/motion distribution
            %is:
            %   - Matrix of no more than 3 dimensions   -> '3d'
            %   - 6 rows exactly                        -> 'nrows', 6
            %   - All values are real, nonnan & numeric -> 'real','nonnan'
            
            if isempty(val) %Always okay
                val = [];
                return
            end
            validateattributes(val, {'numeric'}, {'3d', 'nrows', 6, ...
                'nonnan', 'real'}, class(obj), tok);
        end
    end
    
end

