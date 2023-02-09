classdef FlightPoint < matlab.mixin.SetGet
    %FlightPoint Data class describing the atmospheric conditions at a
    %certain altitude & Mach number that the aircraft is oeprating at.
    %
    % The user can define any combination of flight point data as there is
    % no correlation between any of the atmospheric properties. That is, we
    % do not assume any relation between the gas properties such as the
    % Ideal Gas Law. 
    %
    % The exception to this is if the user only supplies the Mach number
    % and altitude of the flight point. In this instance, when the user
    % requests one of the other air properties (e.g. density, temperature,
    % etc.) it will be calculated using the International Standard
    % Atmosphere equations.
    %
    % The user can invoke the public method 'getFlightPointData' to
    % calculate all of the flight point properties and assign them to the
    % object. Note, this will override any previous data in the object and
    % requries the altitude and Mach number to be defined first. 
    %
    % The altitude must be input in units of FT.
    
    %Aircraft Data
    properties
        %Mach number
        Mach = nan;
        %True air speed (TAS) at the chosen altitude
        AcVelocity
    end
    
    %Flight point data
    properties
        %Altitude [ft]
        Altitude = nan;
        %Air temperature at the chosen altitude
        AirTemp
        %Air density at the chosen altitude
        AirDensity
        %Air pressure at the chosen altitude
        AirPressure
        %Dynamic pressure at the chose altitude (0.5 * rho * V_TAS^2)
        DynPressure
        %Speed of sound at the chosen velocity
        SonicSpeed
        %Ratio of air density at chosen altitude to the reference data
        DensityRatio
    end    
        
    %Sea level values for ISA model
    properties (Constant, Transient)
        %Reference temperature at sea level [K]
        RefTemp     = 288.15;
        %Reference density at sea level [kg/m^3]
        RefDensity  = 1.225;
        %Reference pressure at sea level [N/m^2]
        RefPressure = 101325;
        %Acceleration due to gravity at sea level [m/s^2]
        RefAcc      = 9.80665;
        %Specific gas constant for dry air [J/kgK]
        GasConstant = 287.058;
        %Ratio of specific heats for dry air [-]
        AirGamma    = 1.403;
    end
    
    %List of atmospheric models
    properties (Constant, Transient)
        %Default atmospheric model
        DefaultAtmosphericModel = 'ISA';
        %List of available atmospheric models
        AtmosphericModels       = {'ISA'};
    end
    
    %RefSonicSpeed, Sigma
    properties (Dependent)
        %Speed-of-sound at sea level [m/s]
        RefSonicSpeed
        %Density Ratio
        Sigma
    end
    
    methods % set / get
        function set.Mach(obj, val)           %set.Mach
            %set.Mach Set method for the property 'Mach'
            %
            % Set an arbitrary upper limit of 10 for the Mach number
            
            validateattributes(val, {'numeric'}, {'positive', 'scalar', ...
                'nonempty', 'finite', 'real', '<', 10}, class(obj), 'Mach');
            obj.Mach = val;
        end
        function set.AcVelocity(obj, val)     %set.AcVelocity
            %set.AcVelocity Set method for the property 'AcVelocity'
            
            validateattributes(val, {'numeric'}, {'positive', 'scalar', ...
                'nonempty', 'finite', 'real'}, class(obj), 'AcVelocity');
            obj.AcVelocity = val;
        end
        function set.Altitude(obj, val)       %set.Altitude
            %set.Altitude Set method for the property 'Altitude'
            
            validateattributes(val, {'numeric'}, {'nonnegative', 'scalar', ...
                'nonempty', 'finite', 'real'}, class(obj), 'Altitude'); 
            obj.Altitude = val;
        end
        function set.AirTemp(obj, val)        %set.AirTemp
            %set.AirTemp Set method for the property 'AirTemp'
            
            validateattributes(val, {'numeric'}, {'positive', 'scalar', ...
                'nonempty', 'finite', 'real'}, class(obj), 'AirTemp');
            obj.AirTemp = val;
        end
        function set.AirDensity(obj, val)     %set.AirDensity
            %set.AirDensity Set method for the property 'AirDensity'
            
            validateattributes(val, {'numeric'}, {'positive', 'scalar', ...
                'nonempty', 'finite', 'real'}, class(obj), 'AirDensity');
            obj.AirDensity = val;
        end
        function set.AirPressure(obj, val)    %set.AirPressure
            %set.AirPressure Set method for the property 'AirPressure'
            
            validateattributes(val, {'numeric'}, {'positive', 'scalar', ...
                'nonempty', 'finite', 'real'}, class(obj), 'AirPressure');
            obj.AirPressure = val;
        end
        function set.DynPressure(obj, val)    %set.DynPressure
            %set.DynPressure Set method for the property 'DynPressure'
            
            validateattributes(val, {'numeric'}, {'positive', 'scalar', ...
                'nonempty', 'finite', 'real'}, class(obj), 'DynPressure');
            obj.DynPressure = val;
        end
        function set.SonicSpeed(obj, val)     %set.SonicSpeed
            %set.SonicSpeed Set method for the property 'SonicSpeed'
            
            validateattributes(val, {'numeric'}, {'positive', 'scalar', ...
                'nonempty', 'finite', 'real'}, class(obj), 'SonicSpeed');
            obj.SonicSpeed = val;
        end
    end
    
    methods % calculating the flight point data
        function varargout = getFlightPointData(obj, method, varargin)
            %getFlightPointData Returns the atmospheric conditions at the
            %chosen altitude using the one of the prescribed atmospheric
            %models.

            if nargin < 2 || isempty(method)
                method = obj.DefaultAtmosphericModel;
            end
            
            switch method
                case 'ISA'
                    obj = calculateISA(obj, varargin{:});
                otherwise
                    atm = strjoin(obj(1).AtmosphericModels, ', ');
                    error(['The second input to the method '           , ...
                        '''getFlightPointData'' must be the name of a ', ...
                        'valid atmospheric model.\n\nThe following is ', ...
                        'a list of valid atmospheric model names:'     , ...
                        '\n\n\t%s\n\n'], atm);
            end
            
            if nargout > 0
                varargout{1} = obj;
            end
            
        end        
    end
    
    methods (Access = private) % Atmospheric models
        function obj = calculateISA(obj, varargin)
            %calculateISA Calculates the atmospheric conditions according
            %to the International Standard Atmosphere model.
            
            prp = {'AirTemp', 'AirDensity', 'AirPressure', 'DynPressure', ...
                'AcVelocity', 'SonicSpeed', 'DensityRatio'};
            
            %Parse inputs
            if nargin > 1
                %Identify prop names in varargin
                prp_idx = ismember(prp, varargin);
            else
                prp_idx = true(size(prp));
            end   
            
            if ~any(prp_idx) %Escape route
                return
            end
                        
            ft2m = 0.3048;
            
            %Gather values from (assumed) object arrays
            h      = [obj.Altitude];
            rho0   = [obj.RefDensity];
            T0     = [obj.RefTemp];
            P0     = [obj.RefPressure];
            g0     = [obj.RefAcc];
            R0     = [obj.GasConstant];
            gamma0 = [obj.AirGamma];
            M      = [obj.Mach];
            
            %Check that the altitude has been defined
            h   = h(~isnan(h));
            if isempty(h)
                return
            end
                        
            %Convert the height to m
            h = h .* ft2m;
            
            %There is a change in the variation of temperature and pressure
            %at 11km altitude
            idx = h <= 11000;
            
            %Calculate temperature and pressure for the
            %TROPOSPHERE/STRATOSPHERE region.
            dT = (0.0065 .* h);
            T(idx)  = T0(idx) - dT(idx) ;
            P(idx)  = P0(idx) .* (1 - (dT(idx) ./ T0(idx))) .^ 5.2561;
            T(~idx) = 273 - 56.5;
            P_11    = P0(~idx) .* (1 - (71.5 ./ T0(~idx))) .^ 5.2561;
            P(~idx) = P_11 .* exp( ...
                (-g0(~idx)/(T0(~idx) .* R0(~idx))) .* (h(~idx) - 11000));
            
            %Calculate density, speed of sound, velocity, dynamic pressure
            rho   = P ./ (T .*  R0);
            sigma = rho ./ rho0;
            a     = sqrt(R0 .* gamma0 .* T);
            TAS   = M .* a;
            q     = 0.5 .* rho .* TAS.^2;
            
            %Cell notation
            data{1} = T;
            data{2} = rho;
            data{3} = P;
            data{4} = q;
            data{5} = TAS;
            data{6} = a;
            data{7} = sigma;
            
            %Only assign the information the user has asked for
            prp  = prp(prp_idx);
            data = data(prp_idx);
            arrayfun(@(i) set(obj, prp(i), num2cell(data{i}')), 1 : numel(prp));
            
        end        
    end
    
end

