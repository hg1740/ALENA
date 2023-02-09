classdef MassCase < matlab.mixin.SetGet
    %MassCase Defines a unique combination of payload, fuel and passenger
    %mass and the compartments which the mass will be allocated to.
    
    %Identifying the MassCase
    properties
        %Name of the MassCase
        Name
    end
    
    %Point mass objects specific to this mass case
    properties
        %A vector of 'awi.model.PointMass' objects associated with this
        %mass case
        PointMasses
    end
    
    %Payload mass
    properties
        %Total payload mass (excluding passenger baggage) 
        TotalPayloadMass = 0;
        %Allocation of mass across the different cargo holds
        CargoMasses
        %Allocation of fill-fraction across the different cargo holds
        CargoFraction
        %Handle to the 'awi.model.CargoHold' objects
        CargoHolds
        %Index defining the order which the cargo holds will be loaded
        CargoIndex
    end
    
    %Fuel mass
    properties
        %Total mass of the fuel in all 'awi.model.FuelTank' objects.
        TotalFuelMass = 0;
        %Allocation of fuel across the different fuel tanks
        FuelMasses
        %Allocation of fill-fraction across the different tanks
        FuelFraction
        %Handle to the 'awi.model.FuelTank' objects
        FuelTanks
        %Index defining the order which the fuel tanks will be loaded
        FuelIndex
    end
    
    %Passenger mass
    properties
        %Total mass of all the passengers across the different cabins
        TotalPassengers = 0;
        %Allocation of passengers across the different cabins
        PassengerNumbers
        %Handle to the 'awi.model.Cabin' objects
        Cabins
        %Index defining the order which the cabins will be loaded
        CabinIndex
        %Mass associated with a single passenger
        MassPerPassenger = 0;
        %Baggage mass associated with a single passenger
        CargoMassPerPassenger = 0;
    end    
    
    %Total cargo mass
    properties (Dependent)
        %Total mass in cargo hold (payload + passenger baggage)
        TotalMassInCargoHold
        %Total mass from fuel, cargo and passengers
        TotalMass
    end
    
    methods % set / get
        function set.Name(obj, val)                  %set.Name
           %set.Name Set method for the property 'Name'. 

            validateattributes(val, {'char'}, {}, class(obj), 'Name');
            obj.Name = val;
            
        end
        function set.PointMasses(obj, val)           %set.PointMasses
           %set.PointMasses Set method for the property 'PointMasses'.
           
           if isempty(val)
               obj.PointMasses = [];
               return
           end
           validateattributes(val, {'awi.model.PointMass'}, {'column'}, ...
               class(obj), 'PointMasses');
           obj.PointMasses = val;
           
        end
        function set.TotalPayloadMass(obj, val)      %set.TotalPayloadMass 
            %set.TotalPayloadMass Set method for the property
            %'TotalPayloadMass'.
            
            validateattributes(val, {'numeric'}, {'scalar', 'finite',  ...
                'nonnegative', 'real'}, class(obj), 'TotalPayloadMass');
            obj.TotalPayloadMass = val;
            
        end
        function set.CargoMasses(obj, val)           %set.CargoMasses
            %set.CargoMasses Set method for the property 'CargoMasses'.
            
            validateMasses(obj, val, 'CargoMasses');
            obj.CargoMasses = val;
            
        end
        function set.CargoFraction(obj, val)         %set.CargoFraction 
            %set.CargoFraction Set method for the property
            %'CargoFraction'.
            
            validateFraction(obj, val, 'CargoFraction');
            obj.CargoFraction = val;
            
        end
        function set.CargoHolds(obj, val)            %set.CargoHolds 
            %set.CargoHolds Set method for the property 'CargoHolds'.
            
            validateCompartment(obj, val, 'awi.model.CargoHold', 'CargoHolds');
            obj.CargoHolds = val;
            
        end
        function set.CargoIndex(obj, val)            %set.CargoIndex 
            %set.CargoIndex Set method for the property 'CargoIndex'.
            
            validateIndex(obj, val, 'CargoIndex');
            obj.CargoIndex = val;
            
        end
        function set.TotalFuelMass(obj, val)         %set.TotalFuelMass 
            %set.TotalFuelMass Set method for the property 'TotalFuelMass'.
            
            validateattributes(val, {'numeric'}, {'scalar', 'finite',  ...
                'nonnegative', 'real'}, class(obj), 'TotalFuelMass');
            obj.TotalFuelMass = val;
            
        end
        function set.FuelMasses(obj, val)            %set.FuelMasses
            %set.FuelMasses Set method for the property 'FuelMasses'.
            
            validateMasses(obj, val, 'FuelMasses');
            obj.FuelMasses = val;
            
        end
        function set.FuelFraction(obj, val)          %set.FuelFraction 
            %set.FuelFraction Set method for the property 'FuelFraction'.
            
            validateFraction(obj, val, 'FuelFraction');
            obj.FuelFraction = val;
            
        end
        function set.FuelTanks(obj, val)             %set.FuelTanks
            %set.FuelTanks Set method for the property 'FuelTanks'.
            
            validateCompartment(obj, val, 'awi.model.FuelTank', 'FuelTanks');
            obj.FuelTanks = val;
            
        end
        function set.FuelIndex(obj, val)             %set.FuelIndex 
            %set.FuelIndex Set method for the property 'FuelIndex'.
            
            validateIndex(obj, val, 'FuelIndex');
            obj.FuelIndex = val;
            
        end
        function set.TotalPassengers(obj, val)       %set.TotalPassengers
            %set.TotalPassengers Set method for the property
            %'TotalPassengers'.
            
            validateattributes(val, {'numeric'}, {'scalar', 'integer', ...
                'nonnegative'}, class(obj), 'TotalPassengers');
            obj.TotalPassengers = val;
            
        end
        function set.PassengerNumbers(obj, val)      %set.PassengerNumbers
            %set.PassengerNumbers Set method for the property
            %'PassengerNumbers'.
            
            validateattributes(val, {'numeric'}, {'row', 'integer', ...
                'nonnegative'}, class(obj), 'PassengerNumbers');
            obj.PassengerNumbers = val;
            
        end
        function set.Cabins(obj, val)                %set.Cabins
           %set.Cabins Set method for the property 'Cabins'.
           
           validateCompartment(obj, val, 'awi.model.Cabin', 'Cabins');
           obj.Cabins = val;
           
        end
        function set.CabinIndex(obj, val)            %set.CabinIndex 
            %set.CabinIndex Set method for the property 'FuelIndex'.
            
            validateIndex(obj, val, 'CabinIndex');
            obj.CabinIndex = val;
            
        end
        function set.MassPerPassenger(obj, val)      %set.MassPerPassenger
            %set.MassPerPassenger Set method for the property
            %'MassPerPassenger'.
            
            validateattributes(val, {'numeric'}, {'scalar', 'nonnegative', ...
                'finite', 'real'}, class(obj), 'MassPerPassenger');
            obj.MassPerPassenger = val;
            
        end
        function set.CargoMassPerPassenger(obj, val) %set.CargoMassPerPassenger
            %set.CargoMassPerPassenger Set method for the property
            %'CargoMassPerPassenger'.
           
            validateattributes(val, {'numeric'}, {'scalar', 'nonnegative', ...
                'finite', 'real'}, class(obj), 'CargoMassPerPassenger');
            obj.CargoMassPerPassenger = val;
            
        end
        function val = get.TotalMassInCargoHold(obj) %get.TotalMassInCargoHold
            %get.TotalMassInCargoHold Get method for the dependent property
            %'TotalMassInCargoHold'.
            
            val = obj.TotalPayloadMass + ...
                (obj.TotalPassengers * obj.CargoMassPerPassenger);
            
        end
        function val = get.TotalMass(obj)            %get.TotalMass
           %get.TotalMass Get method for the dependent property 'TotalMass'
           
           cargoMass = obj.TotalMassInCargoHold;
           fuelMass  = obj.TotalFuelMass;
           paxMass   = obj.TotalPassengers * obj.MassPerPassenger;
           val       = cargoMass + fuelMass + paxMass;
           
        end
    end
    
    methods % assigning mass data 
        function allocateMassData(obj, Model, bForce)
            %allocateMassData Assigns the mass data to the various
            %compartments in the model.
            %
            % For the mass data to be useful each compartment must be
            % assigned a fill-fraction (payload & fuel) and each cabin must
            % be assigned a number of passengers.
            
            %Parse
            validateattributes(Model, {'awi.model.Component'}, ...
                {'scalar'}, class(obj), 'Model');
            if nargin < 3
                bForce = false;
            end
            
            %What needs doing?
            bCargo = all(all(~cellfun(@isempty, ...
                get(obj, {'CargoFraction' , 'CargoHolds'}))));
            bFuel  = all(all(~cellfun(@isempty, ...
                get(obj, {'FuelFraction'    , 'FuelTanks'}))));
            bCabin = all(all(~cellfun(@isempty, ...
                get(obj, {'PassengerNumbers', 'Cabins'}))));
            
            if all([bCargo, bFuel, bCabin]) && ~bForce %Escape route
                return
            end
            
            %Find all 'awi.model.Beam' objects 
            %   - Cannot use 'findall' as that only searches for 'Type'
            %     which is overridden by subclasses of 'Beam'.
            flatModel = flatlist(Model);
            Beams = flatModel(arrayfun(@(o) isa(o, 'awi.model.Beam'), flatModel));
            clear flatModel            
            if isempty(Beams) %Escape route
                return
            end
            
            %Collate all 'awi.model.Compartment' objects 
            Compartments = horzcat(Beams.Compartments);
            
            if ~bCargo %Payload    
                i_selectCompartments(obj, Compartments, 'CargoHolds', 'awi.model.CargoHold');
                i_allocatePayloadFraction(obj, 'CargoHolds', ...
                    'TotalMassInCargoHold', 'CargoIndex', ...
                    'CargoMasses', 'CargoFraction');
            end
            if ~bFuel  %Fuel       
                i_selectCompartments(obj, Compartments, 'FuelTanks', 'awi.model.FuelTank');
                i_allocatePayloadFraction(obj, 'FuelTanks', ...
                    'TotalFuelMass', 'FuelIndex', 'FuelMasses', 'FuelFraction');
            end
            if ~bCabin %Passengers 
                i_selectCompartments(obj, Compartments, 'Cabins', 'awi.model.Cabin');
                for ii = 1 : numel(obj)
                    
                    %Determine capacity
                    totalPax       = obj(ii).TotalPassengers;
                    seatsPerCabin  = [obj(ii).Cabins.NumSeats];
                    availableSeats = sum(seatsPerCabin);
                    nCabin         = numel(obj(ii).Cabins);
                    
                    %Check if there is enough total capacity
                    if totalPax > availableSeats
                        error(['Unable to allocate a passenger number ', ...
                            'which is greater than the total number '  , ...
                            'of seats available. Review the requested ', ...
                            'passenger number or adjust the seat layout.']);
                    end
                    
                    %Split equally amongst available compartments?
                    equalPax = totalPax ./ nCabin;
                    
                    %Check if equal split is valid
                    if any(equalPax > seatsPerCabin)
                        %Fill each cabin to capacity until we have
                        %reached the desired passenger numbers
                        if isempty(obj(ii).CabinIndex)
                            obj(ii).CabinIndex = 1 : nCabin;
                        end
                        pax      = seatsPerCabin(obj(ii).CabinIndex);
                        cumSeats = cumsum(pax); 
                        %Find the point where the capacity is reached
                        ind = find(cumSeats >= availableSeats);
                        %Assign remaining mass
                        pax(ind) = availableSeats - cumSeats(ind - 1);
                        pax(ind + 1 : end) = 0;
                    else
                        %Each compartment has an equal mass
                        pax = repmat(equalPax, [1, nCabin]);
                        set(obj(ii), 'CabinIndex', '-equal');
                    end
                    set(obj(ii), 'PassengerNumbers', pax);
                    
                end
                
            end
            
            function i_selectCompartments(obj, Compartments, nam, cls)
                %i_selectCompartments Assigns the handle to the Compartment
                %objects if it has not already been defined.
                %
                % Any MassCases which do not have the relevant Compartment
                % objects assigned will have a reference to all available
                % Compartments of the corresponding type.
                
                Comp = get(obj, {nam});
                idx  = cellfun(@isempty, Comp);
                
                if any(idx)
                    Comp = Compartments(arrayfun(@(o) ...
                        isa(o, cls), Compartments));
                    set(obj(idx), nam, Comp);
                end
                
            end
                        
            function i_allocatePayloadFraction(obj, p_n, p_t, p_i, p_m, p_f)
                %i_allocatePayloadFraction Assigns the fill-fraction for
                %each Compartment.
                %
                %   - 'p_n' : Prop name containing the 'Compartment' handle
                %   - 'p_t' : Prop name of the total mass
                %   - 'p_i' : Prop name of the index variable
                %   - 'p_m' : Prop name of the mass variable
                %   - 'p_f' : Prop name of the fill-fraction variable
                
                %Loop through mass cases and determine amount of payload in
                %each Compartment
                for iO = 1 : numel(obj)
                    
                    %What is the available capacity?
                    cargoVol       = [obj(iO).(p_n).Volume];
                    cargoMassDens  = [obj(iO).(p_n).MassDensity];
                    capacity       = cargoVol .* cargoMassDens;
                    nComp          = numel(capacity);
                    totalCapacity  = sum(capacity);
                    totalCargoMass = obj(iO).(p_t);
                    
                    %Check if there is enough total capacity
                    if totalCargoMass > totalCapacity
                        error(['Unable to allocate a cargo mass which ' , ...
                            'is greater than the total available '      , ...
                            'volume of the cargo hold(s). Review the '  , ...
                            'requested payload mass or adjust the '     , ...
                            'payload density of the cargo hold(s).']);
                    end
                    
                    %Split equally amongst available compartments?
                    equalMass = totalCargoMass ./ nComp;
                    
                    %Check if equal split is valid
                    if any(equalMass > capacity)
                        %Fill each compartment to capacity until we have
                        %reached the desired mass
                        if isempty(obj(iO).(p_i))
                            obj(iO).(p_i) = 1 : nComp;
                        end
                        pmass       = capacity(obj(iO).(p_i));
                        cumCapacity = cumsum(pmass); 
                        %Find the point where the capacity is reached
                        index = find(cumCapacity >= totalCargoMass);
                        %Assign remaining mass
                        pmass(index) = totalCargoMass - cumCapacity(index - 1);
                        pmass(index + 1 : end) = 0;
                    else
                        %Each compartment has an equal mass
                        pmass = repmat(equalMass, [1, nComp]);
                        set(obj(iO), p_i, '-equal');
                    end
                    set(obj(iO), p_m, pmass);
                    
                    %Calculate fill-fraction based on Compartment volume,
                    %mass density and the mass that needs to be allocated
                    fillFrac = arrayfun(@(i) calculatePayloadFraction( ...
                        obj(iO).(p_n)(i), obj(iO).(p_m)(i)), ...
                        1 : numel(obj(iO).CargoHolds), 'Unif', false);
                    
                    %Assign
                    set(obj(iO), p_f, cell2mat(fillFrac));
                    
                end
                
            end
            
        end
    end
    
    methods % converting to FE model
        function MassFEM = convertThisToFE(obj, StructuralFEM)
            %convertThisToFE Converts the 'awi.model.MassCase' object to a
            %collection of Finite Element (FE) entities.
            
            MassFEM = [];
            
            if nargin < 2 %Escape route
                return
            end
            validateattributes(StructuralFEM, {'awi.fe.FEModel'}, ...
                {'nonempty'}, class(obj), 'StructuralFEM');
            
            if numel(obj) > 1 %Object arrays 
                MassFEM = arrayfun(@(o) convertThisToFE(o, StructuralFEM), obj);
                return
            end
                        
            MassFEM      = awi.fe.MassModel;
            MassFEM.Name = obj.Name;
            
            %Generate lumped masses for cargo/fuel compartments
            i_generateFEPointMass(obj.CargoHolds, obj.CargoFraction, MassFEM, StructuralFEM);
            i_generateFEPointMass(obj.FuelTanks , obj.FuelFraction , MassFEM, StructuralFEM);
            
            function i_generateFEPointMass(Comp, pf, MassFEM, StructuralFEM)
                %i_generateFEPointMass Creates the 'awi.fe.PointMass'
                %objects representing the payload in the Compartments
                %'Comp' with payload fraction 'pf'.
                
                %Only proceeed if there is a payload fraction
                Comp = Comp(pf ~= 0);
                
                if numel(Comp) ~= numel(pf) || isempty(Comp) || isempty(pf) %Escape route
                    return
                end
                
                %Handle to geometry objects
                GObj = [StructuralFEM.GeometryObject];
                
                %Define the payload masses
                set(Comp, {'PayloadFraction'}, num2cell(pf)');
                updateMassAndInertia(Comp);
                %Define the 'awi.fe.PointMass' objects for each compartment
                for ii = 1 : numel(Comp)
                    MassProps = Comp(ii).PayloadRigidBodyMassPropsBySection;
                    nMass     = numel(MassProps);
                    FEMass    = arrayfun(@(~) awi.fe.PointMass, 1 : nMass);
                    %Assign data
                    set(FEMass, {'M'}, num2cell([MassProps.Mass])');
                    set(FEMass, {'I11'}, arrayfun(@(mp) ...
                        mp.InertiaTensor(1, 1), MassProps, 'Unif', false)');
                    set(FEMass, {'I21'}, arrayfun(@(mp) ...
                        mp.InertiaTensor(2, 1), MassProps, 'Unif', false)');
                    set(FEMass, {'I22'}, arrayfun(@(mp) ...
                        mp.InertiaTensor(2, 2), MassProps, 'Unif', false)');
                    set(FEMass, {'I31'}, arrayfun(@(mp) ...
                        mp.InertiaTensor(3, 1), MassProps, 'Unif', false)');
                    set(FEMass, {'I32'}, arrayfun(@(mp) ...
                        mp.InertiaTensor(3, 2), MassProps, 'Unif', false)');
                    set(FEMass, {'I33'}, arrayfun(@(mp) ...
                        mp.InertiaTensor(3, 3), MassProps, 'Unif', false)');
                    %Associate 'awi.fe.PointMass' object with a 'Node' on
                    %the parent 'Beam'.
                    FEM = StructuralFEM(ismember(GObj, Comp(ii).BeamHandle));
                    if isempty(FEM) %Escape route
                        continue
                    end
                    assignMassNodes(FEM, FEMass, vertcat(MassProps.CoG));
                    %Add data to the MassFEM
                    addFEData(MassFEM, FEMass);
                end
                
            end
            
            %Generate lumped masses for the passengers  
            i_generateFEPassengerMass(obj.Cabins, obj.PassengerNumbers, ...
                obj.MassPerPassenger, MassFEM, StructuralFEM);
            
            function i_generateFEPassengerMass(Cabin, pax, paxMass, MassFEM, StructuralFEM)
                %i_generateFEPassengerMass Creates the 'awi.fe.PointMass'
                %objects representing the passengers in the Cabins
                %'Cabin' with passenger numbers 'pax' and mass 'paxMass'.
                
                Cabin = Cabin(pax ~= 0);
                
                if numel(Cabin) ~= numel(pax) || isempty(Cabin) || isempty(pax) %Escape route
                    return
                end
                
                %Handle to geometry objects
                GObj = [StructuralFEM.GeometryObject];
                
                %Loop through cabins
                for iC = 1 : numel(Cabin)
                    %How many passengers?
                    nPax  = pax(iC);
                    %Grab the coordinates of the seats
                    xyzS  = Cabin(iC).SeatCoords;
                    %How many seats?
                    nSeat = size(xyzS, 1);
                    %Work out where passengers are sitting
                    if nPax < nSeat
                        %Assume passengers will fill up from the front
                        %   - TODO : We could have all sorts of complex loading
                        %   schemes here!
                        seatsPerRow = Cabin(iC).SeatsPerRow;
                        loadOrder   = cumsum(seatsPerRow);
                        %Find the last full row
                        fullInd    = find(nPax > loadOrder, 1, 'last');
                        if isempty(fullInd)
                            fullInd = 0;
                        end
                        partialInd = fullInd + 1;
                        %Get seat coordinates
                        SeatPosByRow = Cabin(iC).SeatPosByRowInBeamFrame;
                        xFull    = SeatPosByRow.X(1 : fullInd, :);
                        yFull    = SeatPosByRow.Y(1 : fullInd, :);
                        zFull    = SeatPosByRow.Z(1 : fullInd, :);
                        xPartial = SeatPosByRow.X(partialInd, :);
                        yPartial = SeatPosByRow.Y(partialInd, :);
                        zPartial = SeatPosByRow.Z(partialInd, :);
                        %How many seats on the partial row?
                        nParRowSeat = loadOrder(partialInd) - loadOrder(fullInd);
                        nPartial    = nPax - loadOrder(fullInd);
                        nEmpty      = nParRowSeat - nPartial;
                        %Assume the partial row fills up from outside in (i.e.
                        %the window seats fill first then fill inwards)
                        nP = floor(nEmpty / 2);
                        p0 = floor(nParRowSeat / 2);
                        nS = ceil(nEmpty / 2);
                        s0 = ceil(nParRowSeat/2 + 1);
                        pInd = p0 : -1 : p0 - (nP - 1);
                        sInd = s0 :  1 : s0 - (nS - 1);
                        emptyInd = [pInd, sInd];                        
                        %Remove empty seats
                        xPartial(emptyInd) = [];
                        yPartial(emptyInd) = [];
                        zPartial(emptyInd) = [];
                        %Combine with full row coordinates
                        x = [xFull(:) ; xPartial(:)];
                        y = [yFull(:) ; yPartial(:)];
                        z = [zFull(:) ; zPartial(:)];
                        %Gather
                        xyzM = [x, y, z];
                    else
                        %Every seat is filled
                        xyzM = xyzS;
                    end
                    
                    %Convert to global frame
                    xyzM  = xyzM + Cabin(iC).BeamHandle.AbsPosition;
                    
                    %Make the objects
                    PaxMass = arrayfun(@(~) awi.fe.PointMass, 1 : nPax);
                    
                    %Assign data
                    set(PaxMass, 'M', paxMass);
                    
                    %Assign nodes
                    FEM = StructuralFEM(ismember(GObj, Cabin(iC).BeamHandle));
                    if isempty(FEM) %Escape route
                        continue
                    end
                    assignMassNodes(FEM, PaxMass, xyzM);
                    
                    %Add data to the MassFEM
                    addFEData(MassFEM, PaxMass);
                    
                end
                
            end

            %Generate a 'MassSet' object to define the Nastran mass group
%             MassSet = awi.fe.MassSet;
%             MassSet.MassModel = MassFEM;
            addFEData(MassFEM); %, MassSet);            
            
        end
    end
    
    methods (Access = private) % validation
        function validateCompartment(obj, val, cls, prp)
            %validateCompartment Validation method for the generic
            %compartment properties.
            
            validateattributes(val, {cls}, {'row'}, class(obj), prp)
            
        end
        function validateMasses(obj, val, prp)
            %validateMasses Validation method for the generic mass
            %properties.
            
            validateattributes(val, {'numeric'}, {'row', 'finite', ...
                'nonnegative', 'real'}, class(obj), prp);
            
        end
        function validateFraction(obj, val, prp)
            %validateFraction Validation method for the generic
            %fill-fraction properties.
            
            validateattributes(val, {'numeric'}, {'row', 'nonnegative', ...
                '<=', 1, 'finite', 'real'}, class(obj), prp);
            
        end
        function validateIndex(obj, val, prp)
            %validateIndex Validation method for the generic index
            %properties.
            
            if ischar(val)
                validatestring(val, {'-equal'}, class(obj), prp);
            else
                validateattributes(val, {'numeric'}, {'row', 'integer', ...
                    'positive'}, class(obj), prp);
            end
            
        end        
    end
    
end

