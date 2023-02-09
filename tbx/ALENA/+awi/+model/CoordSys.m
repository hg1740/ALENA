classdef (ConstructOnLoad) CoordSys < awi.model.Component & awi.mixin.BeamPropertyObject
    %CoordSys Defines a 3D rotation matrix
    %
    % The rotation matrix can be defined one of 5 ways:
    %   1. Intrinsic Tait-Byran angles
    %   2. Extrinsic Tait-Byran angles
    %   3. Projected sweep, dihedral and AoA (aircraft application only)
    %   4. Quaternions
    %   5. Grid Points
    %   6. Vectors
    %
    % Defining the rotation matrix using one of the above options will
    % cause the other parameters to be back-calculated directly from the
    % rotation matrix.
    %
    % TO DO:
    %   - Add calculation of the rotation vector and angle
    %   - Move the projected set out of the CoordSys as it is too specific
    %   in its application. It should be in the LiftingSurface instead!
    
    %Rotation Matrix (rSet)
    properties (AbortSet, SetObservable)
        %Initial coordinate system before the rotations
        InitialRMatrix = eye(3);
        %3x3 rotation matrix
        RMatrix
    end
    
    %Hidden rotation matrix property
    properties (Dependent, Hidden = true)
        %Rotation matrix in column format - RMatrix_ = RMatrix(:))
        RMatrix_
    end
    
    %Angle rotations (iSet, eSet, pSet)
    properties (AbortSet, SetObservable)
        %Rotation angles about the local axis system (e.g. intrinsic/iSet)
        AngleI
        %Rotation angles about the global axis system (e.g. extrinsic/eSet)
        AngleE
        %Angle between the final rotation matrix axis vectors and the
        %global axes (projected/pSet) - (dih, aoa, swp)
        AngleP
        %Unit of measure for the input angle
        AngleUnits = 'deg';
        %Order of rotations for the 'iSet' and 'eSet'
        RotOrder   = '123';
    end
    
    %Quaternions     (qSet)
    properties (AbortSet, SetObservable)
        %Quaternions
        Quaternions
    end
    
    %Grid set        (gSet)
    properties (AbortSet, SetObservable)
        %A point defining the first vector of the rotation matrix (x,y or
        %z) in the global (x,y,z) system.
        G1
        %A point defining the second vector of the rotation matrix which
        %lies in a plane orthogonal to the first vector.
        G2
        %The vector of the 3x3 rotation matrix which is given by G1
        GridVector = 'z';
        %The plane of the 3x3 rotation matrix which is given by G2
        GridPlane  = 'xz';
    end
    
    %Vector set
    properties (AbortSet, SetObservable)
        %A vector defining the first axis of the rotation matrix (x,y or z)
        V1
        %A vector which is not colinear with 'V1'. The cross product of
        %'V1' and 'V2' defines the 'GridPlane'.
        V2
    end
    
    %Appearance of the coordinate system
    properties (AbortSet, SetObservable)
        ArrowSize
    end
    
    %Valid property values
    properties (Constant, Hidden = true)
        % valid angle units
        ValidAngleUnits    = { ...
            'deg', 'degs', 'degree', 'degrees'; ... % degree set - first character must be 'd'
            'rad', 'rads', 'radian', 'radians'};    % radian set - first character must be 'r'
        % valid rotation order for 'eSet' & 'iSet'
        ValidRotOrder = {'123', '132', '213', '231', '321', '312'}
        % valid grid vectors for the 'gSet'
        ValidGridVector = {'x', 'y', 'z'};
        % valid grid planes for the 'gSet'
        ValidGridPlane = {'xz', 'xy', 'yz'};
    end
    
    %Analytic properties
    properties (Constant, Hidden = true)
        TinyThreshold = 1e-10;   % threshold for determining if a number is near zero
    end
    
    methods % constructor
        function obj = CoordSys(varargin)
            %CoordSys Constructor for the 'CoordSys' class
            
            %Pass it on
            obj@awi.model.Component(varargin{:});
            
            %Add parameter sets
            obj.addParameterSet('rSet', ...
                'DisplayName', 'Rotation Matrix Set', ...
                'Description', 'Define the 3x3 rotation matrix directly.', ...
                'Precedence' , 1, ...
                'RMatrix'    , '3x3 rotation matrix');
            obj.addParameterSet('iSet', ...
                'DisplayName', 'Intrinsic Set', ...
                'Description', 'Intrinsic rotation angles.', ...
                'Precedence' , 2, ...
                'AngleI'     , 'Intrinsic rotation angles about the local x-axis', ...
                'AngleUnits' , 'Units of measure for the angles', ...
                'RotOrder'   , 'Order of rotation');
            obj.addParameterSet('eSet', ...
                'DisplayName', 'Extrinsic Set', ...
                'Description', 'Extrinsic rotation angles.', ...
                'Precedence' , 3, ...
                'AngleE'     , 'Extrinsic rotation about the global x-axis', ...
                'AngleUnits' , 'Units of measure for the angles', ...
                'RotOrder'   , 'Order of rotation');
            obj.addParameterSet('pSet', ...
                'DisplayName', 'Projected Set', ...
                'Description', ['Projected angles between the final axes ', ...
                'of the rotation matrix and the global axes.'], ...
                'Precedence' , 4, ...
                'AngleP'     , 'Projected rotation angle about the global axes', ...
                'AngleUnits' , 'Units of measure for the angles', ...
                'RotOrder'   , 'Order of rotation');
            obj.addParameterSet('qSet', ...
                'DisplayName', 'Quaternion Set', ...
                'Description', 'Quaternions.'  , ...
                'Precedence' , 5, ...
                'Quaternions', 'Quaternion set');
            obj.addParameterSet('gSet', ...
                'DisplayName', 'Grid Set', ...
                'Description', 'Uses coordinates to define a vector and a plane.', ...
                'Precedence' , 6, ...
                'G1'         , '', ...
                'G2'         , '', ...
                'GridVector' , '', ...
                'GridPlane'  , '');
            obj.addParameterSet('vSet', ...
                'DisplayName', 'Vector Set', ...
                'Description', 'Uses two vectors and a cross product to define the rotation matrix', ...
                'Precedence' , 7, ...
                'V1'         , '', ...
                'V2'         , '', ...
                'GridVector' , '', ...
                'GridPlane'  , '');
            
            %CoordSys has no children
            obj.IsLeafNode = true;
            
        end
    end
    
    methods % set / get
        function set.RMatrix(obj, val)     %set.RMatrix     
            %set.RMatrix Set method for the property 'RMatrix'.
            %
            %    - 'RMatrix' must be a numeric of size [3,3] and must be
            %    nonnan, finite and real.
            
            validateattributes(val, {'numeric'}, {'nrows', 3, 'ncols', 3, ...
                'finite', 'real', 'nonnan'}, class(obj), 'RMatrix');
            %normalise column vectors
            obj.RMatrix = cell2mat(arrayfun(@(x)val(:,x)/norm(val(:,x)),...
                1:3,'UniformOutput',false));
            
            %Check for orthoginality
            orthogonalityCheck(obj)
            
            
            
            
        end
        function set.AngleI(obj, val)      %set.AngleI      
            %set.AngleI Set method for the property 'AngleI'
            %
            %    - 'AngleI' must be a numeric row vector with 3 elements and
            %    must be nonnan, finite and real.
            
            validateattributes(val, {'numeric'}, {'row', 'numel', 3, ...
                'nonnan', 'finite', 'real'}, class(obj), 'AngleI')
            obj.AngleI = val;
        end
        function set.AngleE(obj, val)      %set.AngleE      
            %set.AngleE Set method for the property 'AngleE'
            %
            %    - 'AngleE' must be a numeric row vector with 3 elements and
            %    must be nonnan, finite and real.
            
            validateattributes(val, {'numeric'}, {'row', 'numel', 3, ...
                'nonnan', 'finite', 'real'}, class(obj), 'AngleE')
            obj.AngleE = val;
        end
        function set.AngleP(obj, val)      %set.AngleP      
            %set.AngleP Set method for the property 'AngleP'
            %
            %    - 'AngleP' must be a numeric row vector with 3 elements and
            %    must be nonnan, finite and real.
            
            validateattributes(val, {'numeric'}, {'row', 'numel', 3, ...
                'nonnan', 'finite', 'real'}, class(obj), 'AngleP')
            obj.AngleP = val;
        end
        function set.Quaternions(obj, val) %set.Quaternions 
            %set.Quaternions Set method for the property 'Quaternions'
            %
            %    - 'Quaternions' must be a numeric row vector with 4
            %    elements and must be nonnan, finite and real.
            
            validateattributes(val, {'numeric'}, {'row', 'numel', 4, ...
                'nonnan', 'finite'}, class(obj), 'Quaternions'); %, 'real'
            obj.Quaternions = val;
        end
        function set.G1(obj, val)          %set.G1          
            %set.G1 Set method for the property 'G1'.
            %
            %   - 'G1' must be a numeric row vector with 3 elements and
            %    must be nonnan, finite and real.
            
            validateattributes(val, {'numeric'}, {'row', 'numel', 3, ...
                'nonnan', 'finite', 'real'}, class(obj), 'G1')
            obj.G1 = val;
        end
        function set.G2(obj, val)          %set.G2          
            %set.G2 Set method for the property 'G2'.
            %
            %   - 'G2' must be a numeric row vector with 3 elements and
            %    must be nonnan, finite and real.
            
            validateattributes(val, {'numeric'}, {'row', 'numel', 3, ...
                'nonnan', 'finite', 'real'}, class(obj), 'G2')
            obj.G2 = val;
        end
        function set.GridVector(obj, val)  %set.GridVector  
            %set.GridVector Set method for the property 'GridVector'
            %
            % Rules:
            %    - 'GridVector' must be one of the valid strings defined by
            %     'obj.ValidGridVector'
            
            % validate
            validatestring(val, obj.ValidGridVector, class(obj), 'RotOrder');
            % assign --> ensure lower case
            obj.GridVector = lower(val);
        end
        function set.GridPlane(obj, val)   %set.GridPlane   
            %set.GridPlane Set method for the property 'GridPlane'
            %
            % Rules:
            %    - 'GridPlane' must be one of the valid strings defined by
            %     'obj.ValidGridPlane'
            
            % validate
            validatestring(val, obj.ValidGridPlane, class(obj), 'RotOrder');
            % assign --> ensure upper case
            obj.GridPlane = upper(val);
        end
        function set.V1(obj, val)          %set.V1          
            %set.V1 Set method for the property 'V1'
            %
            % Rules:
            %    - 'V1' must be a [1, 3] matrix of real numbers.
            
            validateattributes(val, {'numeric'}, {'row', 'numel', 3, ...
                'real', 'nonnan'}, class(obj), 'V1');
            obj.V1 = val;
        end
        function set.V2(obj, val)          %set.V2          
            %set.V1 Set method for the property 'V2'
            %
            % Rules:
            %    - 'V2' must be a [1, 3] matrix of real numbers.
            
            validateattributes(val, {'numeric'}, {'row', 'numel', 3, ...
                'real', 'nonnan'}, class(obj), 'V2');
            obj.V2 = val;
        end
        function set.AngleUnits(obj, val)  %set.AngleUnits  
            %set.AngleUnits Set method for the property 'AngleUnits'
            %
            % Rules:
            %    - 'AngleUnits' must be one of the valid strings defined by
            %     'obj.ValidAngleUnits'
            
            angleUnits = {'deg', 'rad'};
            
            %Does the token match one of the valid values?
            idx = ismember(obj.ValidAngleUnits', val);
            idx = any(idx);
            
            %Check it matches
            if ~any(idx)
                validatestring(val, obj.ValidAngleUnits(:), class(obj), 'AngleUnits')
            else
                obj.AngleUnits =  angleUnits{idx};
            end
        end
        function set.RotOrder(obj, val)    %set.RotOrder    
            %set.RotOrder Set method for the property 'RotOrder'
            %
            % Rules:
            %    - 'RotOrder' must be one of the valid strings defined by
            %      'obj.ValidRotOrder' BUT if numeric value supplied, value  
            %      is cast to string before validation.
            
            if isnumeric(val)m %Check for numeric
                val = num2str(val);
            end            
            % validate
            validatestring(val, obj.ValidRotOrder, class(obj), 'RotOrder');
            % assign --> ensure upper case
            obj.RotOrder = upper(val);
        end
        function set.ArrowSize(obj, val)   %set.ArrowSize   
            %set.ArrowSize Set method for the property 'ArrowSize'.
            %
            %    - 'ArrowSize' must be a nonnegative, scalar of type
            %    numeric.
            
            validateattributes(val {'numeric'}, {'scalar', 'nonnegative', ...
                'nonnan', 'finite', 'real'}, class(obj), 'ArrowSize');
            obj.ArrowSize = val;
        end
        function val = get.RMatrix(obj)    %get.RMatrix     
            %get.RMatrix Get method for the dependent property 'RMatrix'.
            %
            %    - Default to eye(3) if 'RMatrix' is empty.
            
            if isempty(obj.RMatrix)
                val = eye(3);
            else
                val = obj.RMatrix;
            end
        end
        function val = get.RMatrix_(obj)   %get.RMatrix_    
            %get.RMatrix_ Get method for the property 'RMatrix_'.
            %
            % 'RMatrix_' is the column vector format of the 3x3 rotation
            % matrix.
            
            val = obj.RMatrix(:);
        end
    end
    
    methods % class building
        
        function build_rSet(~)    %build_rSet
            %build_rSet Builds the CoordSys object using the variables
            %belonging to the 'rSet'.
            
            %No action necessary in this method as the 'RMatrix' property
            %is automatically defined by the 'rSet'.
        end
        
        function build_iSet(obj)  %build_iSet
            %build_iSet Builds the CoordSys object using the variables
            %belonging to the 'iSet'.
            
            %Determine basic rotation matricies
            Rx = calcRx(obj, obj.AngleI(1));
            Ry = calcRy(obj, obj.AngleI(2));
            Rz = calcRz(obj, obj.AngleI(3));
            
            %Calculate the rotation matrix based on the order of rotations
            switch obj.RotOrder
                case '123'
                    obj.RMatrix = Rx*Ry*Rz;
                case '132'
                    obj.RMatrix = Rx*Rz*Ry;
                case '213'
                    obj.RMatrix = Ry*Rx*Rz;
                case '231'
                    obj.RMatrix = Ry*Rz*Rx;
                case '321'
                    obj.RMatrix = Rz*Ry*Rx;
                case '312'
                    obj.RMatrix = Rz*Rx*Ry;
            end
            
            obj.RMatrix = obj.RMatrix * obj.InitialRMatrix;
        end
        
        function build_eSet(obj)  %build_eSet
            %build_eSet Builds the CoordSys object using the variables
            %belonging to the 'eSet'.
            
            % determine basic rotation matricies
            Rx = calcRx(obj, obj.AngleE(1));
            Ry = calcRy(obj, obj.AngleE(2));
            Rz = calcRz(obj, obj.Anglee(3));
            
            % calculate the rotation matrix based on the order of rotations
            switch obj.RotOrder
                case '123'
                    obj.RMatrix = Rz*Ry*Rx;
                case '132'
                    obj.RMatrix = Ry*Rz*Rx;
                case '213'
                    obj.RMatrix = Rz*Rx*Ry;
                case '231'
                    obj.RMatrix = Rx*Rz*Ry;
                case '321'
                    obj.RMatrix = Rx*Ry*Rz;
                case '312'
                    obj.RMatrix = Ry*Rx*Rz;
            end
            
            obj.RMatrix = obj.RMatrix * obj.InitialRMatrix;
            
        end
        
        function build_pSet(obj)  %build_pSet
            %build_pSet Builds the CoordSys object using the variables
            %belonging to the 'pSet'.
            
            %Calculate trig. terms
            [~, ~, tX]  = calcTrigValues(obj, obj.AngleP(1));
            [cY, sY, ~] = calcTrigValues(obj, obj.AngleP(2));
            [cZ, sZ, ~] = calcTrigValues(obj, obj.AngleP(3));
            
            %Calculate eY
            eY = [sZ ; cZ ; tX*cZ./(1+tX.^2*cZ.^2).^0.5];
            
            %Calculate eX
            eX_temp = [cY ; 0 ; -sY]; %start with ex projection in XZ plane
            %ex.ey = 0   ==> ex1*ey1 + ex2*ey2 + ex3*ey3 = 0
            %            ==> ey = (-ex1*ey1 - ex3*ey3)/ey2
            eX_temp(2) = (-eX_temp(1)*eY(1) - eX_temp(3)*eY(3))/eY(2);   % create y component of ex to ensure orthogonality with ey
            eX = eX_temp./norm(eX_temp);                                 % normalise to obtain ex vector
            
            %Calculate eZ
            eZ = cross(eX,eY);
            eZ(3) = abs(eZ(3));
            
            %Construct rotation matrix
            obj.RMatrix = [eX, eY, eZ];
            
            obj.RMatrix = obj.RMatrix * obj.InitialRMatrix;
            
        end
                
        function build_qSet(obj)  %build_qSet
            %build_qSet Builds the CoordSys object using the variables
            %belonging to the 'qSet'.
            
            q1 = obj.Quaternion(1);
            q2 = obj.Quaternion(2);
            q3 = obj.Quaternion(3);
            q4 = obj.Quaternion(4);
            
            obj.RMatrix =  [ ...
                q1^2+q2^2-q3^2-q4^2, 2*(q2*q3-q1*q4)    , 2*(q2*q4+q1*q3);
                2*(q2*q3+q1*q4)    , q1^2-q2^2+q3^2-q4^2, 2*(q3*q4-q1*q2);
                2*(q2*q4-q1*q3)    , 2*(q3*q4+q1*q2)    , q1^2-q2^2-q3^2+q4^2];
            
            obj.RMatrix = obj.RMatrix * obj.InitialRMatrix;
            
        end
        
        function build_gSet(obj)  %build_gSet
            %build_gSet Builds the CoordSys object using the variables
            %belonging to the 'gSet'.
            
            tol = 1e-8;
            
            %Parse inputs
            parseGridVectorAndGridPlane(obj);
            
            %Where is the coordinate system in 3D space?
            origin = obj.AbsPosition;
            
            %Calculate the first vector
            switch lower(obj.GridVector)
                case 'x'
                    eX = obj.G1 - origin;
                    eX = eX ./ norm(eX);
                    v1 = eX;
                case 'y'
                    eY = obj.G1 - origin;
                    eY = eY ./ norm(eY);
                    v1 = eY;
                case 'z'
                    eZ = obj.G1 - origin;
                    eZ = eZ ./ norm(eZ);
                    v1 = eZ;
            end
            
            %Vector that lies in the plane defined by 'obj.GridPlane'
            v2 = obj.G2 - origin;
            v2 = v2 ./ norm(v2);
            
            switch lower(obj.GridPlane)
                case {'xy', 'yx'}
                    eZ = cross(v1, v2);
                    eZ = eZ ./ norm(eZ);
                    %Out-of-plane component is always positive
                    if obj.G2(3) + eZ(3) - origin(3) < tol
                        eZ = eZ - (2 .* eZ);
                    end
                    %Now define the final vector
                    switch lower(obj.GridVector)
                        case 'x'
                            eY = cross(eX, eZ);
%                             eY = checkVectorDirection(eY, v2, 'xy');
%                             if sign(eY(1)) ~= sign(v2(1))
%                                 eY = eY - (2 .* eY);
%                             end
                        case 'y'
                            eX = cross(eY, eZ);
                    end
                case {'xz', 'zx'}
                    eY = cross(v1, v2);
                    eY = eY ./ norm(eY);
                    %Out-of-plane component is always positive
                    if obj.G2(2) + eZ(2) < origin(2)
                        eY = eY - (2 .* eY);
                    end
                    switch lower(obj.GridVector)
                        case 'x'
                            eZ = cross(eX, eY);
                        case 'z'
                            eX = cross(eY, eZ);
                    end
                case {'yz', 'zy'}
                    eX = cross(v1, v2);
                    eX = eX ./ norm(eX);
                    %Out-of-plane component is always positive
                    if obj.G2(1) + eX(1) < origin(1)
                        eX = eX - (2 .* eX);
                    end
                    switch lower(obj.GridVector)
                        case 'y'
                            eZ = cross(eX, eY);
                        case 'z'
                            eY = cross(eX, eZ);
                    end
            end
            
            %Construct the rotation matrix
            obj.RMatrix = [eX', eY', eZ'];
            
            function v = checkVectorDirection(v, v0, plane)
                
                %Which quadrant
                switch plane
                    case 'xy'
                        ind = [1, 2];
                    case 'yz'
                        ind = [2, 3];
                    case 'xz'
                        ind = [1, 3];
                end
                q  = getQuadrant(v, ind);
                q0 = getQuadrant(v0, ind);
                
                if q ~= q0
                    v = mapVector2Quadrant(v, ind, q, q0);
                end
                
            end
            
            function q = getQuadrant(v, ind)
                %getQuadrant Returns which qudrant the vector is in.
                
                v_1 = v(ind(1));
                v_2 = v(ind(2));
                
                if v_1 >= 0 && v_2 >= 0
                    q = '1';
                elseif v_1 < 0 && v_2 >= 0 
                    q = '2';
                elseif v_1 < 0 && v_2 < 0
                    q = '3';
                elseif v_1 >= 0 && v_2 < 0
                    q = '4';
                end
                
            end
            
            function v = mapVector2Quadrant(v, ind, q1, qN)
                %mapVector2Quadrant Maps a single vector from one quadrant
                %to another
                
                if q1 == qN
                    return
                end
                
                if q1 == '1' && qN == '2'
                    v(ind(2)) = abs(v(ind(2)));                    
%                 elseif q1 == '1' &7 qN
                    
                end
                
            end
            
        end
        
        function build_vSet(obj)  %build_vSet
            %build_vSet Builds the CoordSys object using the variables
            %belonging to the 'vSet'.
            
            tol = 1e-8;
            
            %Parse inputs
            parseGridVectorAndGridPlane(obj);            

            %Normalise the vectors
            v1 = obj.V1 ./ norm(obj.V1);
            v2 = obj.V2 ./ norm(obj.V2);
            
            %Calculate the first vector
            switch lower(obj.GridVector)
                case 'x'
                    eX = v1;
                case 'y'
                    eY = v1;
                case 'z'
                    eZ = v1;
            end
            
            %Calculate the out-of-plane vector 
            switch lower(obj.GridPlane)
                case {'xy', 'yx'}
                    eZ = cross(v1, v2);
                    eZ = eZ ./ norm(eZ);
%                     %Out-of-plane component is always positive
%                     if v2(3) + eZ(3) < tol
%                         eZ = eZ - (2 .* eZ);
%                     end
                    %Now define the final vector
                    switch lower(obj.GridVector)
                        case 'x'
                            eY = cross(eX, eZ);
%                             eY = checkVectorDirection(eY, v2, 'xy');
%                             if sign(eY(1)) ~= sign(v2(1))
%                                 eY = eY - (2 .* eY);
%                             end
                        case 'y'
                            eX = cross(eY, eZ);
                    end
                case {'xz', 'zx'}
                    eY = cross(v1, v2);
                    eY = eY ./ norm(eY);
%                     %Out-of-plane component is always positive
%                     if v2(2) + eY(2) < tol
%                         eY = eY - (2 .* eY);
%                     end
                    switch lower(obj.GridVector)
                        case 'x'
                            eZ = cross(eX, eY);
                        case 'z'
                            eX = cross(eY, eZ);
                    end
                case {'yz', 'zy'}
                    eX = cross(v1, v2);
                    eX = eX ./ norm(eX);
%                     %Out-of-plane component is always positive
%                     if v2(1) + eX(1) < tol
%                         eX = eX - (2 .* eX);
%                     end
                    switch lower(obj.GridVector)
                        case 'y'
                            eZ = cross(eX, eY);
                        case 'z'
                            eY = cross(eX, eZ);
                    end
            end
            
            %Construct the rotation matrix
            obj.RMatrix = [eX', eY', eZ'];
            
        end
        
        function update_rSet(~)   %update_rSet
            %update_rSet Updates the 'rSet' based on the current values of
            %'RMatrix'.
            
            %No action necessary as 'RMatrix' is already defined by all
            %other build methods.
            
        end
        
        function update_iSet(obj) %update_iSet
            %update_iSet Updates the 'iSet' based on the current values of
            %'RMatrix'.
            
            switch obj.AngleUnits
                case 'deg'
                    ia(1) = atan2d(obj.RMatrix(3,2), obj.RMatrix(3,3));
                    ia(2) = asind(-obj.RMatrix(3,1));
                    ia(3) = atan2d(obj.RMatrix(2,1), obj.RMatrix(1,1));
                case 'rad'
                    ia(1) = atan2(obj.RMatrix(3, 2), obj.RMatrix(3,3));
                    ia(2) = asin(-obj.RMatrix(3,1));
                    ia(3) = atan2(obj.RMatrix(2,1), obj.RMatrix(1,1));
            end
            
            obj.AngleI = ia;
            
        end
        
        function update_eSet(obj) %update_eSet
            %update_eSet Updates the 'eSet' based on the current values of
            %'RMatrix'.
            
            switch obj.AngleUnits
                case 'deg'
                    ea(1) = atan2d(-obj.RMatrix(2,3), obj.RMatrix(3,3));
                    ea(2) = asind(obj.RMatrix(1,3));
                    ea(3) = atan2d(-obj.RMatrix(1,2), obj.RMatrix(1,1));
                case 'rad'
                    ea(1) = atan2(-obj.RMatrix(2,3), obj.RMatrix(3,3));
                    ea(2) = asin(obj.RMatrix(1,3));
                    ea(3) = atan2(-obj.RMatrix(1,2), obj.RMatrix(1,1));
            end
            
            obj.AngleE = ea;
            
        end
        
        function update_pSet(obj) %update_pSet
            %update_pSet Updates the 'pSet' based on the current values of
            %'RMatrix'.
            
            % define global coordinate axes
            xG = [1,0,0];
            yG = [0,1,0];
            zG = [0,0,1];
            
            % define normal vectors to global coordinate planes
            XY = cross(xG, yG);
            XZ = cross(zG, xG);
            YZ = cross(yG, zG);
            
            % grab local axis vectors from the rotation matrix
            eX = obj.RMatrix(:, 1, :)';
            eY = obj.RMatrix(:, 2, :)';
            eZ = obj.RMatrix(:, 3, :)';
            
            % get vector projections
            eX_XZ = obj.calcVectorProjection(eX, XZ);
            eY_XY = obj.calcVectorProjection(eY, XY);
            eY_YZ = obj.calcVectorProjection(eZ, YZ);
            
            % angleX is the angle between the global Y-axis and the vector eY
            pa(1) = obj.calcVectorAngle(eY_YZ, XZ, obj.AngleUnits);
            pa(2) = obj.calcVectorAngle(eY_XY, XZ, obj.AngleUnits);
            pa(3) = obj.calcVectorAngle(eX_XZ, YZ, obj.AngleUnits);
            
            obj.AngleP = pa;
        end
        
        function update_qSet(obj) %update_qSet
            %update_qSet Updates the 'qSet' based on the current values of
            %'RMatrix'.
            
            %             error('Check the Rxx, Rxy, Rxz, etc. terms');
            
            % grab terms from the rotation matrix
            Rxx = obj.RMatrix(1,1); Rxy = obj.RMatrix(1,2); Rxz = obj.RMatrix(1,3);
            Ryx = obj.RMatrix(2,1); Ryy = obj.RMatrix(2,2); Ryz = obj.RMatrix(2,3);
            Rzx = obj.RMatrix(3,1); Rzy = obj.RMatrix(3,2); Rzz = obj.RMatrix(3,3);
            
            w = sqrt(trace(obj.RMatrix) + 1) / 2;
            
            % check if w is real. Otherwise, zero it.
            if( imag(w) > 0 )
                w = 0;
            end
            
            x = sqrt(1 + Rxx - Ryy - Rzz) / 2;
            y = sqrt(1 + Ryy - Rxx - Rzz) / 2;
            z = sqrt(1 + Rzz - Ryy - Rxx) / 2;
            
            [~, i] = max([w, x, y, z]);
            
            if( i == 1 )
                x = ( Rzy - Ryz ) / (4*w);
                y = ( Rxz - Rzx ) / (4*w);
                z = ( Ryx - Rxy ) / (4*w);
            end
            
            if( i == 2 )
                w = ( Rzy - Ryz ) / (4*x);
                y = ( Rxy + Ryx ) / (4*x);
                z = ( Rzx + Rxz ) / (4*x);
            end
            
            if( i == 3 )
                w = ( Rxz - Rzx ) / (4*y);
                x = ( Rxy + Ryx ) / (4*y);
                z = ( Ryz + Rzy ) / (4*y);
            end
            
            if( i == 4 )
                w = ( Ryx - Rxy ) / (4*z);
                x = ( Rzx + Rxz ) / (4*z);
                y = ( Ryz + Rzy ) / (4*z);
            end
            
            % assign quaternion to the object
            q(1) = w;
            q(2) = x;
            q(3) = y;
            q(4) = z;
            
            obj.Quaternions = q;
            
        end
        
        function update_gSet(obj) %update_gSet
            %update_gSet Updates the 'gSet' based on the current values of
            %'RMatrix'.
            
            %First grid point is one of the rotation matrix axes
            switch obj.GridVector
                case 'x'
                    obj.G1 = obj.RMatrix(:, 1)';
                case 'y'
                    obj.G1 = obj.RMatrix(:, 2)';
                case 'z'
                    obj.G1 = obj.RMatrix(:, 3)';
            end
            
            %Second grid point lies in the same plane as the second point
            switch obj.GridPlane
                case 'xz'
                    if strcmp(obj.GridVector, 'x')
                        obj.G2 = obj.RMatrix(:, 3)';
                    elseif strcmp(obj.GridVector, 'z')
                        obj.G2 = obj.RMatrix(:, 1)';
                    end
                case 'xy'
                    if strcmp(obj.GridVector, 'x')
                        obj.G2 = obj.RMatrix(:, 2)';
                    elseif strcmp(obj.GridVector, 'y')
                        obj.G2 = obj.RMatrix(:, 1)';
                    end
                case 'yz'
                    if strcmp(obj.GridVector, 'y')
                        obj.G2 = obj.RMatrix(:, 3)';
                    elseif strcmp(obj.GridVector, 'z')
                        obj.G2 = obj.RMatrix(:, 2)';
                    end
            end
            
            %Make sure we are in the global system!
            obj.G1 = obj.G1 + obj.AbsPosition;
            obj.G2 = obj.G2 + obj.AbsPosition;
            
        end
        
        function update_vSet(obj) %update_vSet 
            
%             error('Update code!!!');
        end
        
    end
    
    methods % visualisation
        
        function hg = drawElement(obj, ht, tag)
            %drawElement Class specific draw method for the 'CoordSys'
            %class.
            
            %Caller supplied tag explicitly ?
            if nargin < 3
                
                %No - apply a tag indicating that we are just drawing some origins
                tag = 'Coordinate Systems';
                
            end
                               
            %Grab principal axes of rotation matrix
            R  = cat(3, obj.RMatrix);
            eX = squeeze(R(:, 1, :));
            eY = squeeze(R(:, 2, :));
            eZ = squeeze(R(:, 3, :));
            
            %Grab origin and 
            O  = vertcat(obj.Position)';
            eX = eX + O;
            eY = eY + O;
            eZ = eZ + O;
            
            %Axis vectors - Pad with NaN to allow vectorisation and
            %minimise number of hg objects
%             Ox   = O(:, 1); Oy = O(:, 2); Oz = O(:, 3);
%             nan_ = nan(size(Ox));            
%             OXx = [Ox, eX(:, 1), nan_];
%             OXy = [Oy, eX(:, 2), nan_];
%             OXz = [Oz, eX(:, 3), nan_];
%             OYx = [Ox, eY(:, 1), nan_];
%             OYy = [Oy, eY(:, 2), nan_];
%             OYz = [Oz, eY(:, 3), nan_];
%             OZx = [Ox, eZ(:, 1), nan_];
%             OZy = [Oy, eZ(:, 2), nan_];
%             OZz = [Oz, eZ(:, 3), nan_];
            Ox   = O(1, :); Oy = O(2, :); Oz = O(3, :);
            nan_ = nan(size(Ox));            
            OXx = [Ox ; eX(1, :) ; nan_];
            OXy = [Oy ; eX(2, :) ; nan_];
            OXz = [Oz ; eX(3, :) ; nan_];
            OYx = [Ox ; eY(1, :) ; nan_];
            OYy = [Oy ; eY(2, :) ; nan_];
            OYz = [Oz ; eY(3, :) ; nan_];
            OZx = [Ox ; eZ(1, :) ; nan_];
            OZy = [Oy ; eZ(2, :) ; nan_];
            OZz = [Oz ; eZ(3, :) ; nan_];

            %Set in column format
            OXx = OXx(:); OXy = OXy(:); OXz = OXz(:);
            OYx = OYx(:); OYy = OYy(:); OYz = OYz(:);
            OZx = OZx(:); OZy = OZy(:); OZz = OZz(:);
            
            %X-axis in red
            hg{1} = line('Parent', ht, ...
                'XData', OXx, ...
                'YData', OXy, ...
                'ZData', OXz, ...
                'Marker', 'none', ...
                'LineStyle', '-', ...
                'Color', 'r', ...
                'Tag', tag);
            
            %Y-axis in green
            hg{2} = line('Parent', ht, ...
                'XData', OYx, ...
                'YData', OYy, ...
                'ZData', OYz, ...
                'Marker', 'none', ...
                'LineStyle', '-', ...
                'Color', 'g', ...
                'Tag', tag);
            
            %Z-axis in blue
            hg{3} = line('Parent', ht, ...
                'XData', OZx, ...
                'YData', OZy, ...
                'ZData', OZz, ...
                'Marker', 'none', ...
                'LineStyle', '-', ...
                'Color', 'b', ...
                'Tag', tag);
            
            %Grab coordinates at the end of each axis vector
            OXTxt = [OXx(2 : 3 : end), OXy(2 : 3 : end), OXz(2 : 3 : end)];
            OYTxt = [OYx(2 : 3 : end), OYy(2 : 3 : end), OYz(2 : 3 : end)];
            OZTxt = [OZx(2 : 3 : end), OZy(2 : 3 : end), OZz(2 : 3 : end)];
            
            %Construct coordinates, label & color for the Text object
            txtCoords = [OXTxt ; OYTxt ; OZTxt];
            txtString = repmat({'X' , 'Y' , 'Z'}, [numel(obj), 1]);
            txtString = txtString(:);
            txtColor  = repmat({'r', 'g', 'b'}, [numel(obj), 1]);
            txtColor  = txtColor(:);
            
            %Add the text
            hT = text(txtCoords(:, 1), txtCoords(:, 2), txtCoords(:, 3), ...
                txtString, ...
                'Parent'    , ht    , ...
                'Tag'       , tag   , ...
                'FontWeight', 'bold', ...
                'SelectionHighlight', 'off');
            
            set(hT, {'Color'}, txtColor);
%             set(hT, 'FontWeight', 'bold');
%             set(hT, 'SelectionHighlight', 'off');
            set([hg{:}], 'SelectionHighlight', 'off');
            
        end
        
    end
    
    methods % trigonometric functions 

        function [c, s, t] = calcTrigValues(obj, theta)
            %calcTrigValues Calculates the sin, cos and tan of some angle
            %'theta' accounting for whether the angle is provided in degrees or
            %radians
            
            switch obj.AngleUnits
                case 'deg'
                    c = cosd(theta);
                    s = sind(theta);
                    t = tand(theta);
                case 'rad'
                    c = cos(theta);
                    s = sin(theta);
                    t = tan(theta);
            end
        end

        function Rx = calcRx(obj, theta)
            %calcRx Calculates the basic 3D rotation matrix for a clockwise
            %rotation about the x-axis by some angle 'theta'
            
            [c, s, ~] = calcTrigValues(obj, theta);
            Rx = [1, 0 ,0 ; 0, c, -s ; 0, s, c ];
        end

        function Ry = calcRy(obj, theta)
            %calcRy Calculates the basic 3D rotation matrix for a clockwise
            %rotation about the y-axis by some angle 'theta'
            
            [c, s, ~] = calcTrigValues(obj, theta);
            Ry = [c, 0, s ; 0, 1, 0 ; -s, 0, c];
        end

        function Rz = calcRz(obj, theta)
            %calcRz Calculates the basic 3D rotation matrix for a clockwise
            %rotation about the z-axis by some angle 'theta'
            
            [c, s, ~] = calcTrigValues(obj, theta);
            Rz = [c, -s, 0 ; s, c, 0 ; 0, 0, 1];
        end
        
    end
    
    methods (Static) % calcVectorAngle, calcVectorProjection, skew, rodriguezRotation

        function angle = calcVectorAngle(u, v, angleUnits)
            %calcVectorAngle Calculates the angle between two vectors.
            %
            % It is assumed that the vectors 'u' and 'v' are defined such that
            % the columns denote the values in the global (x,y,z) coordinate
            % system.
            %
            % e.g. u = [ux, uy, uz], v = [vx, vy, vz]
            %
            % This function can handle a vectorised input as long as
            
            % parse inputs
            validateattributes(u, {'numeric'}, {'nonempty', 'nonnan', 'finite', 'real', 'ncols', 3, '2d'});
            validateattributes(v, {'numeric'}, {'nonempty', 'nonnan', 'finite', 'real', 'ncols', 3, '2d'});
            validatestring(angleUnits, {'deg', 'rad'}, 'calcVectorAngle', 'angleUnits');
            
            % check 'u' & 'v' are the same size
            if ~isequal(size(u), size(v))
                ME = MException('MATLAB:uob:CoordSys:sizemismatch', ...
                    'The size of the vectors ''u'' and ''v'' must be the same');
                throw(ME);
            end
            
            % calculate angle between the two vectors
            switch angleUnits
                case 'deg'
                    angle = atan2d(norm(cross(u,v,2)),dot(u,v,2));
                case 'rad'
                    angle = atan2(norm(cross(u,v,2)),dot(u,v,2));
            end
        end

        function uP = calcVectorProjection(u, v)
            %calcVectorProjection Calculates the projection of the vector 'u'
            %onto the plane or vector defined by 'v'.
            %
            %More information can be found <a
            % href="http://www.euclideanspace.com/maths/geometry/elements/plane/lineOnPlane/">here</a>.
            
            uP = cross(v,cross(u,(v./norm(v,2))./norm(v,2),2),2);
        end
        
        function skewMatrix = skew(vector)
            %skew Returns the skew-symmetric matrix of a vector.
            
            skewMatrix   = [ ...
                0         ,-vector(3),  vector(2) ; ...
                 vector(3),    0     , -vector(1) ; ...
                -vector(2), vector(1),    0      ];
            
        end
        
        function rotMatrix = rodriguezRotation(rot_vector, rot_angle)
            %rodriguezRotation Returns the rotation matrix resulting from a
            %rotation of 'rot_angle' degrees about the vector 'rot_vector'.
            %
            %<a href="https://en.wikipedia.org/wiki/Rodrigues%27_rotation_formula">Wikipedia page for Rodriguez Rotation Formula</a>
            
            K    = awi.model.CoordSys.skew(rot_vector);
            sina = sind(rot_angle);
            cosa = cosd(rot_angle);
            rotMatrix = eye(3) + (sina * K) + (1 - cosa) * K^2;
            
        end
        
    end
    
    methods (Access = private) % helper functions
        function parseGridVectorAndGridPlane(obj)
            %parseGridVectorAndGridPlane Parses the values of 'GridVector'
            %and 'GridPlane' to ensure that the combination of the two will
            %create a valid orthogonal rotation matrix. 
            
            switch lower(obj.GridPlane)
                case {'xy', 'yx'}
                    if strcmpi(obj.GridVector, 'z')
                        error(['Cannot define the coordinate system for ', ...
                            'the case where G1 defines the z-axis and '  , ...
                            'G2 defines a point in the XY-plane']);
                    end
                case {'xz', 'zx'}
                    if strcmpi(obj.GridVector, 'y')
                        error(['Cannot define the coordinate system for ', ...
                            'the case where G1 defines the y-axis and '  , ...
                            'G2 defines a point in the XZ-plane']);
                    end
                case {'yz', 'zy'}
                    if strcmpi(obj.GridVector, 'x')
                        error(['Cannot define the coordinate system for ', ...
                            'the case where G1 defines the x-axis and '  , ...
                            'G2 defines a point in the YZ-plane']);
                    end
            end
            
        end
    end
    
    % error handling
    methods 
        function orthogonalityCheck(obj)
            %orthogonalityCheck Checks for orthogonality between the
            %(x,y,z) vectors in the rotation matrix.
            
            %Grab axis vectors
            eX = obj.RMatrix(:, 1, :);
            eY = obj.RMatrix(:, 2, :);
            eZ = obj.RMatrix(:, 3, :);
            
            %Check for othogonality between vectors of the rotation matrix
            if (dot(eX,eY) > obj.TinyThreshold) %eX & eY
                ME = obj.generateOrthME('eX', 'eY', obj);
                throwAsCaller(ME);
            end
            if (dot(eX,eZ) > obj.TinyThreshold) %eX & eZ
                ME = obj.generateOrthME('eX', 'eZ', obj);
                throwAsCaller(ME);
            end
            if (dot(eY,eZ) > obj.TinyThreshold) %eY & eZ
                ME = obj.generateOrthME('eY', 'eZ', obj);
                throwAsCaller(ME);
            end
        end
    end    
    methods (Static) 
        function ME = generateSetME(setName, setProps)
            %generateSetME Generates the MException object for handling
            %incorrect inputs for the constructor method of the 'CoordSys'
            %class.
            
            propString = strjoin(setProps, ', ');
            nArg       = length(setProps);
            % define MException object
            ME = MException('MATLAB:uob:CoordSys:wrongInputs', ...
                ['When defining the rotation matrix using the ''%s'' '   , ...
                'option you must provide %i additional arguments. These ', ...
                'are %s. For more information, see the %s documention.'] , ...
                setName, nArg, propString, 'CoordSys');
            % throw error to user
            throwAsCaller(ME);
        end        
        function ME = generateOrthME(vec1, vec2, obj)
            %generateOrthME Generates the MExeption object for when the axis
            %vectors of the rotation matrix are not orthogonal.
            
            ME = MException('MATLAB:uob:CoordSys:vectorsNotNormal', ...
                ['The %s and %s vectors of the calculated rotation '    , ...
                'matrix are not orthogonal.\n\nThe dot product of the ' , ...
                'vectors %s and %s should be less than '                , ...
                '''%s.TinyThreshold'' (%.8g).'], vec1, vec2, vec1, vec2 , ...
                class(obj), obj.TinyThreshold);
        end        
    end
    
end


