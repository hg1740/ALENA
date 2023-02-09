classdef (ConstructOnLoad) LoadCase < awi.model.Entity
    %LoadCase represents loading on an aero-elastic body
    
    properties (SetObservable) 
        %Identifier
        ID = NaN;        
        %Mass of the payload
        PayloadMass     = 0;
        %Fraction of the available fuel volume that is being used
        FuelFraction    = 0;
        %Mass of the fuel associated with this loadcase
        FuelMass        = 0;
        %Number of passengers
        NumPassenger    = 0;
        %Load factor in the vertical (positive lift) direction
        LoadFactor      = 1;
        %Safety factor
        SafetyFactor    = 1.5;
        %Centre-of-Gravity as a percentage of the wing MAC
        CgMac           = 0.25;
        %Aircraft pitch angle
        PitchAngle      = 0;
        %Aircraft roll angle
        RollAngle       = 0;
        %Thrust of a single engine
        SingleEngineThrust = 0;
        %Mach number
        Mach            = NaN;
        %Altitude
        Altitude        = NaN;
        %True air speed (TAS) at the chosen altitude
        AcVelocity      = NaN;
        %Aircraft All-Up-Mass (AUM)
        AcMass          = 0;
        %Control surface(s) have a deflection, and a "type" (e.g. fixed, or free)
        CsDeflection;
        CsDeflecType = {};        
        %For convenience, control surface deflection properties can be tabulated
        CsDeflections;        
        %To allow distinction between manouevre, 1MC and CT load cases
        LoadCaseType;
        %GRAV card flag
        GRAV = 9.81;
        %Rigid Trim Flag
        Flexibility = 1;
        %CS Flag
        AileronAngle = 0;
        %CS Flag
        ElevatorAngle = 0;
        %CS Flag
        RudderAngle = 0;
    end
        
    %Gust load case properties
    properties
        %Applicable to 1MC and CT load case types
        GustLength = 10:10:110;
        FrontDirection;
        GustDirection;        
        %Appicable to CT load case type only
        TurbulenceScale;
        SpectrumType;
    end
    
    %Cruise Design Case properties
    properties
        %Applicable to Cruise case only
        TargetCl_eta;
        TargetCl;        
    end
    
    %Helper functions for populating the LoadCase options
    properties (Constant, Transient) 
        %Types of load cases that can be defined
        LoadCaseTypes = {'Manoeuvre', 'One-Minus-Cosine', 'Continuous-Turbulence', 'Cruise Design Case', 'Pratt Gust', 'Static'};
        %Types of atmospheric turbulence spectrums that are accepted
        SpectrumTypes = {'Von Karman', 'Dryden'};
        %Allowable state(s) for the control surfaces
        CsDeflecTypes = {'fixed', 'free'};        
    end
    
    %Control surface parameters
    properties (Dependent)
        %Handle to all 'awi.model.ControlSurface' objects in the Framework
        ControlSurfaces
        %Names of all the control surfaces in the Framework
        ControlSurfaceNames
        %Handle to all of the parents of the individual control surfaces
        ControlSurfaceParents
        %Handle to all of the parent beams of the individual control surfaces
        ControlSurfaceParentBeam
    end
    
    %Flight point data
    properties 
        FlightPoint = []; %TODO - Swap out all air properties and contain in a 'awi.model.FlightPoint' object instead
        %Air temperature at the chosen altitude
        AirTemp
        %Air density at the chosen altitude
        AirDensity
        %Air pressure at the chosen altitude
        AirPressure
        %Dynamic pressure at the choseN altitude (0.5 * rho * V_TAS^2)
        DynPressure
        %Speed of sound at the chosen velocity
        SonicSpeed
    end
    
    %Sea level values for ISA mode
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
    properties (Dependent)
        %Speed-of-sound at sea level [m/s]
        RefSonicSpeed
    end
   
    %Fuel distribution
    properties
        %Handle to a vector of 'awi.model.FuelDistribution' objects
        FuelDistribution
    end
    
    %Mass case data
    properties
        %Handle to an 'awi.model.MassCase' object
        MassCase = awi.model.MassCase.empty;
        %Mass of a single passenger
        PassengerMass
        %Mass of hold luggage
        LuggageMass
    end    
    
    methods % set / get
        function set.LoadFactor(obj, val)
            
            %Cannot set the load factor for a cruise method 
            %   - TODO Phil to set property to locked in the GUI
            assert(val == 1 || ~strcmp(obj.LoadCaseType, 'Cruise Design Case'), ...
                'Load factor cannot be set for a load case of type ''Cruise Design Case''.') %#ok<MCSUP>
            
            obj.LoadFactor = val;
            
        end
        
        function val = get.LoadFactor(obj)
            
            val = obj.LoadFactor;
            
            %A cruise load case must have a load factor of 1.
            if strcmp(obj.LoadCaseType, 'Cruise Design Case')
                val = 1;
            end
        end

        function set.GRAV(obj, val)

            obj.GRAV = val;

        end
        function val = get.GRAV(obj)

            val = obj.GRAV;

        end
        function set.Flexibility(obj, val)

            obj.Flexibility = val;

        end
        function val = get.Flexibility(obj)

            val = obj.Flexibility;

        end
        function val = get.AileronAngle(obj)

            val = obj.AileronAngle;

        end
        function val = get.ElevatorAngle(obj)

            val = obj.ElevatorAngle;

        end
        function val = get.RudderAngle(obj)

            val = obj.RudderAngle;

        end
        function val = get.SpectrumType(obj)
        
            %Start here
            val = obj.SpectrumType;
            
            %Nothing yet ?
            if isempty(val)
                val = obj.SpectrumTypes{1};
            end
            
        end
    
        function set.SpectrumType(obj, val)
            
            %Validate
            assert(ismember(val, obj.SpectrumTypes), ['invalid spectrum type : ', val]); 
            
            %Make a note
            obj.SpectrumType = val;
            
        end
            
        function val = get.LoadCaseType(obj)
        
            %Start here
            val = obj.LoadCaseType;
            
            %Nothing yet ?
            if isempty(val)
                val = obj.LoadCaseTypes{1};
            end
            
        end
    
        function set.LoadCaseType(obj, val)
            
            %Validate
            assert(ismember(val, obj.LoadCaseTypes), ['invalid load case type : ',val]); 
            
            %Make a note
            obj.LoadCaseType = val;
            
        end
        
        function val = get.ControlSurfaces(obj)
            %TODO - This will fail if we ever have more than one aircraft
            %in the Framework hierachy
            if isa(obj.Root, 'awi.model.Framework')
                
                % -> Chris I changed this as it would otherwise ignore
                % control surfaces gathered by the collector
                %val = obj.Root.ControlSurfaces;
                Collection = flatlist(obj.Root);
                cidx   = arrayfun(@(o) isa(o, 'awi.model.ControlSurface'), Collection);
                val    = Collection(cidx);
            else
                val = [];
           end
        end
        
        function val = get.ControlSurfaceNames(obj)
            cs = obj.ControlSurfaces;
            if isempty(cs)
                val = [];
            else
                val = {cs.Name};
            end
        end
        
        function val = get.ControlSurfaceParents(obj)
            cs = obj.ControlSurfaces;            
            if isempty(cs)
                val = [];
            else
                val = [cs.Parent]; %TODO - Parents might not be set
            end
        end
        
        function val = get.ControlSurfaceParentBeam(obj)
            cs = obj.ControlSurfaces;            
            if isempty(cs)
                val = [];
            else
                val = [cs.ParentBeam]; %TODO - Parents might not be set
            end
        end
        
        function val = get.CsDeflections(obj)
            
            %Get control surface names
            csNames   = obj.ControlSurfaceNames;
            
            %Get parents
            csParents = obj.ControlSurfaceParentBeam;
            
            %Get parent names - allowing for parents being not set
            if isempty(csParents)
                csPName = {};
            else
                csPName = {csParents.Name};
            end
            
            %Get current value(s)
            csValues = obj.CsDeflection;
            
            %Pad with zeros
            csValues((end + 1) : (end+(numel(csNames) - numel(csValues)))) = 0;
            
            %Remove surplus data
            csValues(numel(csNames) + 1 : end) = [];
           
            %Ditto for control surface types,but infill with defaults as we go
            csTypes = obj.CsDeflecType;
            csTypes(cellfun(@isempty, csTypes)) = repmat(obj.CsDeflecTypes(1), 1, sum((cellfun(@isempty, csTypes))));
            csTypes((end + 1) : (end+(numel(csNames) - numel(csTypes)))) = repmat(obj.CsDeflecTypes(1), 1, numel(csNames) - numel(csTypes));
            csTypes(numel(csNames) + 1 : end) = [];
            
            %Return as a table
            val = table(csNames(:), csPName(:), csValues(:), csTypes(:), ...
                'VariableNames', {'Name', 'Parent', 'Value', 'Type'});
            
            %Use UserData to carry detail of enumerations around with table (TODO: good idea ?)
            val.Properties.UserData{4} = obj.CsDeflecTypes;
            
        end
        
        function set.FuelDistribution(obj, val)
            %set.FuelDistribution Set method for the property
            %'FuelDistribution'.
            %
            %   - 'FuelDistribution' must be a vector of
            %     'awi.model.FuelDistribution' objects.
            
            %Setting the value to empty
            if isempty(val)
                obj.FuelDistribution = [];
                return
            end
            
            validateattributes(val, {'awi.model.FuelDistribution'}, ...
                {'row', 'nonempty'}, class(obj), 'FuelDistribution');
            obj.FuelDistribution = val;
            
        end
        
        function val = get.AirTemp(obj)       %get.AirTemp       
            %get.AirTemp Get method for the dependent property 'AirTemp'
            %
            % If the air temperature is undefined then use the method
            % 'getFlightPointData' to calculate the air temperature at the
            % specified altitude.
            
            val = obj.AirTemp;
            
            if isempty(val) && ~isnan(obj.Altitude)
               FP  = getFlightPointData(obj);
               val = FP.Temp;
            end
            
        end        
        function val = get.AirDensity(obj)    %get.AirDensity    
            %get.AirDensity Get method for the dependent property 
            %'AirDensity'
            %
            % If the air density is undefined then use the method
            % 'getFlightPointData' to calculate the air density at the
            % specified altitude.   
            
            val = obj.AirDensity;   
            
            if isempty(val) && ~isnan(obj.Altitude)
               FP  = getFlightPointData(obj);
               val = FP.Density;
            end            
            
        end        
        function val = get.AirPressure(obj)   %get.AirPressure   
            %get.AirPressure Get method for the dependent property 
            %'AirPressure'
            %
            % If the air pressure is undefined then use the method
            % 'getFlightPointData' to calculate the air pressure at the
            % specified altitude.   
            
            val = obj.AirPressure;   
            
            if isempty(val) && ~isnan(obj.Altitude)
               FP  = getFlightPointData(obj);
               val = FP.Pressure;
            end            
            
        end
        function val = get.AcVelocity(obj)    %get.AcVelocity    
           %get.AcVelocity Get method for the property 'AcVeclocity' 
           %
           % If the aircraft forward velocity (TAS) is undefined then use
           % the method 'getFlightPointData' to calculate it at the
           % specified altitude.
           
           val = obj.AcVelocity;   
            
            if isempty(val) && ~any(isnan([obj.Altitude, obj.Mach])) 
                FP = getFlightPointData(obj);
               val = FP.Velocity;
            end           
            
        end
        function val = get.DynPressure(obj)   %get.DynPressure   
           %get.DynPressure Get method for the property 'DynPressure' 
           %
           % If the dynamic pressure is undefined then use the method 
           %'getFlightPointData' to calculate it at the specified altitude.
           
           val = obj.DynPressure;   
            
            if isempty(val) && ~any(isnan([obj.Altitude, obj.Mach])) 
               FP = getFlightPointData(obj);
               val = FP.DynPressure;
            end           
            
        end
        function val = get.SonicSpeed(obj)    %get.SonicSpeed    
            %get.SonicSpeed Get method for the property 'SonicSpeed' 
           %
           % If the speed of sound is undefined then use the method 
           %'getFlightPointData' to calculate it at the specified altitude.
           
           val = obj.SonicSpeed;   
            
            if isempty(val) && ~isnan(obj.Altitude)
               FP = getFlightPointData(obj);
               val = FP.SonicSpeed;
            end           
            
        end
        function val = get.RefSonicSpeed(obj) %get.RefSonicSpeed 
           %get.RefSonicSpeed Get method for the dependent property
           %'RefSonicSpeed'.
           
           val = sqrt(obj.GasConstant * obj.AirGamma * obj.RefTemp);
           
        end
        
        function set.MassCase(obj, val)       %set.MassCase 
            %set.MassCase Set method for the property 'MassCase'.
            
            validateattributes(val, {'awi.model.MassCase'}, {'scalar'}, ...
                class(obj), 'MassCase');
            obj.MassCase = val;
            
        end
        function set.NumPassenger(obj, val)   %set.NumPassenger 
            %set.NumPassenger Set method for the property 'NumPassenger'.
            
            validateattributes(val, {'numeric'}, {'scalar', 'integer', ...
                'nonnegative', 'finite', 'real'}, class(obj), ...
                'NumPassenger');
            obj.NumPassenger = val;
            
        end
        function set.PassengerMass(obj, val)  %set.PassengerMass
            %set.PassengerMass Set method for the property 'PassengerMass'.
            
            validateattributes(val, {'numeric'}, {'scalar', 'nonnegative', ...
                'finite', 'real'}, class(obj), 'PassengerMass');
            obj.PassengerMass = val;
            
        end
        function set.LuggageMass(obj, val)    %set.LuggageMass
            %set.LuggageMass Set method for the property 'LuggageMass'.
            
            validateattributes(val, {'numeric'}, {'scalar', 'nonnegative', ...
                'finite', 'real'}, class(obj), 'LuggageMass');
            obj.LuggageMass = val;
            
        end
        
    end
    
    methods % construction / destruction
        
        function obj = LoadCase(varargin)
            
            %TODO - Add units for the different properties.
            
            %Pass it on
            obj@awi.model.Entity(varargin{:});
            
            %Load case has no chldren
            obj.IsLeafNode = true;
            
            %Load cases do not draw themselves in a Drawing view
            obj.Visible = false;
            
            %If properties are managed by Nameable
            if isa(obj, 'mvc.mixin.Nameable')
                
                %Extend property groups
                obj.addPropertyGroup('Properties', ...
                    'LoadCaseType' , 'Load Case Type'    , [], obj.LoadCaseTypes, ...
                    'ID'           , 'Case ID'           , [], [], ...
                    'Altitude'     , 'Altitude'          , [], [], ...
                    'LoadFactor'   , 'Load Factor'       , [], [], ...
                    'FuelFraction' , 'Fuel Fraction'     , [], [], ...
                    'Mach'         , 'Mach'              , [], [], ...
                    'PayloadMass'  , 'Payload Mass'      , [], [], ...
                    'CgMac'        , 'CG position (%MAC)', [], [],...
                    'AcMass'       , 'Aircraft Mass'     , [], [],...
                    'CsDeflections', 'Control Surfaces'  , [], []);
                
                %Properties applicable to 'One-Minus-Cosine'
                obj.addPropertyGroup('One-Minus-Cosine', ...
                    'GustLength'   , 'Gust Gradient(s)', 'Gust Gradient - Distance to peak gust velocity', [], ...
                    'GustDirection', 'Gust Direction'  , 'Gust Direction', []);
                
                %Properties applicable to 'Continuous Turbulence'
                obj.addPropertyGroup('Continuous Turbulence', ...
                    'GustDirection'  , 'Gust Direction'  , 'Gust Direction'   , [], ...
                    'TurbulenceScale', 'Turbulence Scale', 'Turbulence Scale' , [], ...
                    'SpectrumType'   , 'Spectrum Type'   , 'Spectrum type'    , obj.SpectrumTypes);
                
                %Properties applicable to 'Cruise Design Case'
                obj.addPropertyGroup('Cruise Design Case', ...
                    'TargetCl_eta' , 'Target Cl Eta'         , 'Target Cl Eta'         , [], ...
                    'TargetCl'     , 'Target Cl Distribution', 'Target Cl Distribution', []);                
                
                %Configure tab visibilities
                obj.setPropertyGroup('One-Minus-Cosine', ...
                    'Visible', @()strcmp(obj.LoadCaseType, obj.LoadCaseTypes{2}));
                obj.setPropertyGroup('Continuous Turbulence', ...
                    'Visible', @()strcmp(obj.LoadCaseType, obj.LoadCaseTypes{3}));
                obj.setPropertyGroup('Cruise Design Case', ...
                    'Visible', @()strcmp(obj.LoadCaseType, obj.LoadCaseTypes{4}));

                %Configure enable state
                obj.setPropertyAttribute('Properties', 'LoadFactor', ...
                    'Enable', @()strcmp(obj.LoadCaseType, obj.LoadCaseTypes{1}));
                obj.setPropertyAttribute('Properties', 'CsDeflections', ...
                    'Enable', @()ismember(obj.LoadCaseType, obj.LoadCaseTypes([1, 4])));
                
%                 obj.setPropertyAttribute('Properties', 'TargetCl_eta', ...
%                     'Enable', @()~strcmp(obj.LoadCaseType, obj.LoadCaseTypes{1}));
%                 
%                 obj.setPropertyAttribute('Properties', 'TargetCl', ...
%                     'Enable', @()~strcmp(obj.LoadCaseType, obj.LoadCaseTypes{1}));
                 
                %Configure how we make CsDeflections editable
                obj.setPropertyAttribute('Properties', 'CsDeflections', 'Extras', ...
                    {'ColumnEditable', [false, false, true, true], ...
                    'CellEditCallback', @obj.onCsDeflectionsCellEdit});
                
            end
            
        end
        
    end
    
    methods 
       
        function ControlSurf = getControlSurfaceInfo(obj, filter) 
            %getControlSurfaceInfo Retrieves the handles to the control
            %surface objects in the model whose state corresponds to the
            %string in 'filter'.
            
            %Handle reduced inputs
            if nargin < 2
                filter = 'all';
            end
            
            %Validate the filter string
            validatestring(fiter, {'all', 'fixed', 'free'}, class(obj), 'filter');
            
            ControlSurf = obj.ControlSurfaces;
            
            
        end
        
        function FlightPoint = getFlightPointData(obj)
            %getFlightPointData Returns the atmospheric conditions at the
            %load case altitude using the formulas for the International
            %Standard Atmosphere (ISA).
            %   - The altitude is assumed to be in units of 'ft'
            %   - This function supports object arrays.
            %   - The equations for the ISA model are take from:
            %       "The International Standard Atmoshpere (ISA)" - Mustafa
            %       Cavcar, 2000.
            
            FlightPoint = [];
            ft2m = 0.3048;
            
            %Gather values from (assumed) object arrays
            h      = [obj.Altitude];
            T0     = [obj.RefTemp];
            P0     = [obj.RefPressure];
            rho0   = [obj.RefDensity];
            g0     = [obj.RefAcc];
            R0     = [obj.GasConstant];
            gamma0 = [obj.AirGamma];
            M      = [obj.Mach];
            
            %Check that the altitude has been defined            
            h   = h(~isnan(h));
            if isempty(h)
                return
            end
            
            %Preallocate
            prp = {'Temp', 'Density', 'Pressure', ...
                'SonicSpeed', 'AcVelocity', 'DynPressure', 'DensityRatio'}; 
            FlightPoint = repmat(cell2struct(cell(size(prp)), prp, 2), size(obj));
            
            %Convert the height to m
            h = h .* ft2m;
            
            %There is a change in the variation of temperature and pressure
            %at 11km altitude
            idx = h <= 11000;
            
            %Calculate temperature and pressure for the
            %TROPOSPHERE/STRATOSPHERE region.
            dT = (0.0065 .* h(idx));
            T(idx)  = T0(idx) - dT ;
            P(idx)  = P0(idx) .* (1 - (dT ./ T0(idx))) .^ 5.2561;
            T(~idx) = 273 - 56.5;
            P_11    = P0(~idx) .* (1 - (71.5 ./ T0(~idx))) .^ 5.2561;
            P(~idx) = P_11 .* exp( ...
                (-g0(~idx)/(T0(~idx) .* R0(~idx))) .* (h(~idx) - 11000));
            
            %Calculate density, speed of sound, velocity, dynamic pressure
            rho   = P ./ (T .*  R0);
            a     = sqrt(R0 .* gamma0 .* T);
            TAS   = M .* a;
            q     = 0.5 .* rho .* TAS.^2; 
            sigma = rho ./ rho0;
            
            %Cell notation
            T     = num2cell(T);
            P     = num2cell(P);
            rho   = num2cell(rho);
            a     = num2cell(a);
            TAS   = num2cell(TAS);
            q     = num2cell(q);
            sigma = num2cell(sigma);
            M     = num2cell(M);
            
            %Assign to the structure
            [FlightPoint.Temp]         = deal(T{:});
            [FlightPoint.Pressure]     = deal(P{:});
            [FlightPoint.Density]      = deal(rho{:});
            [FlightPoint.DensityRatio] = deal(sigma{:});
            [FlightPoint.SonicSpeed]   = deal(a{:});
            [FlightPoint.AcVelocity]   = deal(TAS{:});
            [FlightPoint.DynPressure]  = deal(q{:});
            [FlightPoint.Mach]         = deal(M{:});
        end
        
        function generateMassCases(obj, Aircraft, massCaseName, bForce)
            %generateMassCases Creates a series of 'awi.model.MassCase'
            %objects based on the mass data in the 'awi.model.LoadCase'
            %objects.
            %
            % For each MassCase object the total payload, fuel and
            % passenger mass is split amongst all available compartments.
            %
            % If a more detailed mass case is desired then the user should
            % create their own and manually assign it to the LoadCase.            
                        
            %Parse
            validateattributes(Aircraft, {'awi.model.Aircraft'}, ...
                {'scalar'}, class(obj), 'generateMassCases');
            if nargin < 3
                massCaseName = {};
            end
            if nargin < 4
                bForce = false;
            end
            
            %Preserve any existing mass cases
            MassCases = {obj.MassCase};
            idx       = cellfun(@isempty, MassCases);
            
            if bForce %Can override if we want to! 
                idx = true(size(idx));
            end
            
            %If the LoadCase objects already have MassCase data then there
            %is no need to proceed unless the user is forcing new mass data
            %to be generated.
            if all(~idx)
                return
            end
            
            %Retrieve mass data
            massPrps = {'PassengerMass', 'LuggageMass', 'PayloadMass', ...
                'FuelMass', 'NumPassenger'};            
            massData = cell2mat(get(obj(idx), massPrps));
            
            if all(massData == 0) %Escape route
                return
            end
            
            %Filter for unique combinations only - Defines no. mass cases
            [uMassData, ia, ib] = unique(massData, 'rows');
            uMassData = num2cell(uMassData);
            nMassCase = numel(ia);
            
            %Name?
            if isempty(massCaseName) || numel(massCaseName) ~= nMassCase
                massCaseName = arrayfun(@(i) sprintf('MassCase%i', i), ...
                    1 : nMassCase, 'Unif', false);
            end
            
            %Define 'awi.model.MassCase' objects and assign data
            MassCases = arrayfun(@(~) awi.model.MassCase, 1 : nMassCase);
            set(MassCases, {'MassPerPassenger'}     , uMassData(:, 1));
            set(MassCases, {'CargoMassPerPassenger'}, uMassData(:, 2));
            set(MassCases, {'TotalPayloadMass'}     , uMassData(:, 3));
            set(MassCases, {'TotalFuelMass'}        , uMassData(:, 4));
            set(MassCases, {'TotalPassengers'}      , uMassData(:, 5));
            set(MassCases, {'Name'}                 , massCaseName');
            
            %Allocate mass data to the various Compartments
            allocateMassData(MassCases, Aircraft, bForce);
            
            %Assign unique MassCases to their respective LoadCases
            MC = cell(size(obj));
            MC = MC(idx);
            for i = 1 : numel(ia)
                MC(ib == i) = {MassCases(i)};
            end
            set(obj(idx), {'MassCase'}, MC);
            
        end
        
        function varargout = calculateGustLoadFactor(obj, Aircraft, CL_alfa)
            %calculateGustLoadFactor Calculates the vertical load factor
            %using the equation for a Pratt Gust.                        
            
            validateattributes(Aircraft, {'awi.model.Aircraft'}, ...
                {'scalar'}, 'calculateGustLoadFactor', 'Aircraft');
            if nargin < 3
                %Assume the aircraft lift-curve slope is equal to 2pi
                CL_alfa = 2 * pi;
            else
                validateattributes(CL_alfa, {'numeric'}, {'row', 'positive'}, 'calculateGustLoadFactor', 'CL_alfa');
            end
            
            %Only proceed using LoadCases of type 'Pratt Gust'
            obj = obj(ismember({obj.LoadCaseType}, 'Pratt Gust'));
            if isempty(obj) %Escape route
                return
            end
            
            %Check we have the correct inputs
            props = get(Aircraft, {'RefArea', 'RefSpan'});
            if any(cellfun(@isempty, props)) %Escape route 
                return
            end
            area = props{1};
            span = props{2};
            chrd = area / span;
            props = get(obj, {'Altitude', 'AcMass', 'Mach'});
            if any(any(cellfun(@isempty, props))) || ...
                    any(any(cellfun(@isnan, props))) %Escape route 
                return
            end
            FP   = getFlightPointData(obj);
            rho  = [FP.Density];
            FPv  = [FP.AcVelocity];
            alt  = horzcat(props{:, 1});
            mass = horzcat(props{:, 2});
            
            %Sometimes the aircraft velocity is explicitly defined. If that
            %is the case then use this data instead of the flight point
            %data.
            vel      = [obj.AcVelocity];
            idx      = isnan(vel);
            vel(idx) = FPv(idx);
            
            %Get mass ratio 'muG' and gust load alleviation factor 'kG'
            muG = mass ./ (0.5 .* rho .* area .* chrd .* CL_alfa);
            kG  = (0.88 .* muG) ./ (5.3 + muG);
            
            %Calculate gust velocity
            Uref = [ ...
                0    , 15000, 60000 ; ...
                17.07, 13.41, 6.36 ];
            Ug = interp1(Uref(1, :), Uref(2, :), alt);  %[m/s], EAS
            Ug = Ug ./ sqrt(rho ./ obj(1).RefDensity);  %[m/s], TAS
            
            %Calculate incremental load factor
            nZ = 1 + (kG .* CL_alfa .* rho .* area .* vel .* Ug ./ (2 .* mass * 9.81));           
            %VT = kG .*(1 + Ug ./ vel) .* (1 + Ug ./ vel) .^2;
            
            %User asked for it back?
            if nargout > 0
                varargout{1} = nZ;
            else
                %Assign to the object
                set(obj, {'LoadFactor'}, num2cell(nZ'));
            end            
            
        end
        
        function acVelocity = getAcVelocity(obj)
            %getAcVelocity Retrieves the aircraft forward velocity.
            
            %Velocity according to flight point
            FP   = getFlightPointData(obj);
            FPv  = [FP.AcVelocity];
            
            %Sometimes the aircraft velocity is explicitly defined. If that
            %is the case then use this data instead of the flight point
            %data.
            acVelocity      = [obj.AcVelocity];
            idx             = isnan(acVelocity);
            acVelocity(idx) = FPv(idx);
            
        end
        
    end
    
    methods (Access = protected) % callbacks for uitable
        
        function onCsDeflectionsCellEdit(obj, ~, evtdata)
            
            %Careful
            try
                
                %Get index of element being edited
                idx = evtdata.Indices(1);
                
                %Which field if being edited ?
                if evtdata.Indices(2) == 3
                    
                    %Deflection itself - get new value
                    val = str2num(evtdata.NewData); %#ok<ST2NM>
                    
                    %Assign in object
                    obj.CsDeflection(idx) = val;
                    
                elseif evtdata.Indices(2) == 4
                    
                    %Deflection type - get new value
                    val = evtdata.NewData;
                    
                    %Assign in object
                    obj.CsDeflecType{idx} = val;
                    
                else
                    error('unhandled column');
                end
                
            catch err
                
                %What went wrong ?
                uiwait(errordlg(obj, err, 'modal'));
                
                %Need to revert value, I guess ? (TODO)
                
            end
            
        end
        
    end
    
    methods (Static) % fame2obj
       
        function obj = fame2obj(lc, lc_fm4, fuelcases)
            %fame2obj Converts a FAME load case structure into an instance
            %of 'awi.model.LoadCase'.
            %
            % Variable description:
            %   - 'lc' : Structure with load case data that has been taken
            %   from the 'load_cases.txt' file in the FAME output.
            %   - 'lc_fm4', : Structure with load case data that has been
            %   taken from the .fm4 input file. This structure will not
            %   have the 1g design cruise case defined.
            %   - 'fuelcases' : Objects descibing the fuel configurations
            %   generated using the files provided by the user.
            %   
            
            %Token/string for matching data in the FAME loadcase structure
            %to the AWI object.
            tokStr = { ...
                'LoadFactor'  , 'Loadfactor'  ; ...
                'SafetyFactor', 'Safetyfactor'; ...
                'CgMac'       , 'AC_cg'       ; ...
                'PitchAngle'  , 'PitchAngle'  ; ...
                'RollAngle'   , 'Rollangle'   ; ...                
                'Mach'        , 'Mach'        ; ...
                'Altitude'    , 'Altitude'    ; ...
                'PayloadMass' ,'PayloadMass'  ; ...
                'SingleEngineThrust' , 'EngineThrust' ; ...
                'AcMass'      , 'AC_mass'};
            
            %Default to no fuelcases provided
            if nargin < 3
                fuelcases = [];
            end
            
            %How many loadcases?
            nLC = numel(lc);
            
            %Create the objects
            obj = arrayfun(@(i) awi.model.LoadCase, 1 : nLC, 'Unif', false);
            obj = horzcat(obj{:});
            
            %Set all those properties that map straight across
            for iTok = 1 : numel(tokStr) / 2
                set(obj, tokStr(iTok, 1), {lc.(tokStr{iTok, 2})}');                
            end
                       
           %Default all load cases to be of type 'Manoeuvre' except those
           %that have a load factor of 1.
           %    - TODO : What if there is more than one 1g case? Better to
           %    just select the end loadcase as it is almost guaranteed to
           %    be the one that is attached by FAME?
           idx_design_case = ismember([lc.Loadfactor], 1);
           set(obj, 'LoadCaseType' , 'Manoeuvre');
           set(obj(idx_design_case), 'LoadCaseType', 'Cruise Design Case');
           set(obj(idx_design_case), 'Name', sprintf('LoadCase %i - Design Cruise Case', numel(obj)));
           set(obj(idx_design_case), 'ID'  , numel(obj));
           
           %Assign data from the .fm4 load cases
           %    - ID           -> ID
           %    - LABEL & TYPE -> Name
           %    - COMMENT      -> Description
           if ~isempty(lc_fm4)
               %Filter the 'awi.model.LoadCase' objects as the .fm4 load
               %case data does not have any data for the 1g design case
               obj_ = obj(~idx_design_case);
               %ID
               id = {lc_fm4.ID};
               id = cellfun(@(x) str2double(x), id, 'Unif', false);
               set(obj_, {'ID'}, id');
               %Name - Combine ID & TYPE (Should give unique)
               type  = {lc_fm4.TYPE};               
               nam   = arrayfun(@(i) sprintf('LoadCase %i - %s', id{i}, type{i}), 1 : numel(type), 'Unif', false);
               set(obj_, {'Name'}, nam');
               %Description - Combine LABEL & COMMENT
               lab  = {lc_fm4.LABEL};
               com  = {lc_fm4.COMMENT};
               %Check for embedded cells
               idxL = cellfun(@(x) iscell(x), lab);
               idxC = cellfun(@(x) iscell(x), com);
               lab(idxL) = cellfun(@(x) x{:}, lab(idxL), 'Unif', false);
               com(idxC) = cellfun(@(x) x{:}, com(idxC), 'Unif', false);
               %Construct the 'Description'
               desc = arrayfun(@(i) [lab{i}, ', ', com{i}], 1 : numel(lab), 'Unif', false);
               set(obj_, {'Description'}, desc');
           end
           
           %Calculate fuel-fraction and assign to loadcase objects
           fuel = [lc.fuel];
           fuelfraction = num2cell(fuel / max(fuel))';
           set(obj, {'FuelFraction'}, fuelfraction);
            
           %Grab name of fuel cases
           if ~isempty(fuelcases)
               %Names
               fuelcaseNames = cellfun(@(x) x(1).FameFuelFile, fuelcases, 'Unif', false);
               %Fuel mass
               %   - Assume the fuel file has been named using the
               %   convention 'SomeStringHere_FuelMass'. This will allow
               %   us to delimit by '_'
               fuelValue = cellfun(@(x) sscanf(x(strfind(x, '_') + 1 : end), '%g'), fuelcaseNames, 'Unif', false);
               %Check everything has been extracted correctly
               idx = cellfun(@(x) isempty(x), fuelValue);
               %TODO - Need to force this to work
               if any(idx)
                   return
               end
               %                 assert(~any(idx), 'Unknown fuel file format');
               %Convert to vector
               fuelValue = horzcat(fuelValue{:});
           end
           
           %If fuel cases have been provided then assign them
           for iFC = 1 : numel(fuelcases)
               %Find the corresponding fuel case based on the
               %mass/volume of the fuel. TODO - Finesse this!
               diff = abs((fuel - fuelValue(iFC)));
               %idx  = diff < 1; %Use a tolerance of 1kg
               idx = (diff == min(diff));
               set(obj(idx), 'FuelDistribution', fuelcases{iFC});
           end           

        end
                
    end
    
end
