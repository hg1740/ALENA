classdef cpacs < awi.model.Entity
    %cpacs Handles the import and export of CPACS models to the AWI
    %Framework.
    %
    % CPACS (Common Parametric Aircraft Schema) is an aircraft
    % parameterisation schema developed by DLR. It is open-source and is
    % accessible via GitHub. For further details see ...
    %
    % The CPACS parameterisation can handle aircraft and rotorcraft,
    % however the AWI framework will only import AIRCRAFT models.
            
    %Function handle for logging progress and informing the user
    properties
        %Function handle to the logging process, by default any log
        %messages will be printed to the command window.
        LogFcn = @(s) fprintf('%s\n', s);
    end
    
    %Flags for controlling the conversion process
    properties (Hidden)
        %Indicates whether the conversion process should he ended.
        EndConversion  = false;
        %Indicates whether the conversion process should be exited.
        QuitConversion = false;
    end
    
    %Handles to the components of the CPACS model
    properties (Dependent)
        %Handle to all of the VEHICLE elements in the CPACS model
        Vehicles
        %Handle to all of the AIRCRAFT elements in the CPACS model
        Aircraft
        %Handle to all of the PROFILE elements in the CPACS model
        Profiles
        %Handle to all of the FUSELAGE PROFILE elements in the CPACS model
        FuselageProfiles
        %Handle to all of the WING AIRFOIL elements in the CPACS model
        WingProfiles
    end
    
    %Shortcuts
    properties (Dependent)
        %Logical flag indicating whether the CPACS model has any vehicles
        HasVehicles
        %Logical flag indicating whether the CPACS model has any aircraft
        %as children of the vehicle.
        HasAircraft
        %Logical flag indicating whether the CPACS model has any profile
        %information as children of the vehicle
        HasProfiles
    end
    
    methods % set / get
        function set.LogFcn(obj, val)            %set.LogFcn           
            validateattributes(val, {'function_handle'}, {'scalar'}, ...
                class(obj), 'LogFcn');
            obj.LogFcn = val;
        end
        function val = get.Vehicles(obj)         %get.Vehicles         
            %get.Vehicles Returns a handle to the vehicle(s) in the CPACS
            %model.
            val = findall(obj, 'Name', 'vehicles');
            if isempty(val)
                val = [];
            end                
        end
        function val = get.Aircraft(obj)         %get.Aircraft         
            %get.Aircraft Returns the aircraft elements within the CPACS
            %model. Could potentially be an array of objects.
            val  = findall(obj.Vehicles, 'Name', 'aircraft');
            if isempty(val)
                val = [];
            end
        end     
        function val = get.Profiles(obj)         %get.Profiles         
            val = findall(obj.Vehicles, 'Name', 'profiles');
            if isempty(val)
                val = [];
            end
        end
        function val = get.FuselageProfiles(obj) %get.FuselageProfiles 
            val = findall(obj.Profiles, 'Name', 'fuselageProfiles');
            if isempty(val)
                val = [];
            end
        end
        function val = get.WingProfiles(obj)     %get.WingProfiles     
            val = findall(obj.Profiles, 'Name', 'wingAirfoils');
            if isempty(val)
                val = [];
            end
        end
        function val = get.HasVehicles(obj)      %get.HasVehicles      
            if isempty(obj.Vehicles)
                val = false;
            else
                val = true;
            end            
        end
        function val = get.HasAircraft(obj)      %get.HasAircraft      
            if isempty(obj.Aircraft)
                val = false;
            else
                val = true;
            end            
        end
        function val = get.HasProfiles(obj)      %get.HasProfiles      
            if isempty(obj.Profiles)
                val = false;
            else
                val = true;
            end     
        end
    end    
        
    methods (Sealed)
        function AWIObj = convertCPACsToAWI(obj, logfcn)
            %convertCPACSToAWI Converts a hierachy of CPACS objects into a
            %valid AWI aircraft model.
            %
            % 'obj' is of class 'awi.model.cpacs' which is a very thin
            % wrapper around 'awi.model.Entity'.            
            
            if nargin < 1 %Make a note of the log function
                obj.LogFcn = logfcn;
            end
            
            %Start with something sensible
            AWIObj = [];
            
            %% Extract the data from the raw CPACS xml output
            
            while ~obj.EndConversion                
                
                obj.LogFcn('Converting CPACS data to AWI Framework model...');
                
                if ~obj.HasVehicles %Check for vehicles 
                    terminateConversion(obj, 'vehicles')
                end
                if ~obj.HasAircraft %Check for aicraft  
                    terminateConversion(obj, 'aircraft');
                end
            
                %How many aicraft do we have?
                ac  = obj.Aircraft;
                nAC = numel(ac);   
                
                %Extract the data for each aircraft (vehicle) in the model
                obj.LogFcn('Extracting the vehicle (aircraft) data...');
                AircraftData = arrayfun(@(iAC) obj.extractAircraftData( ...
                    ac(iAC), obj.LogFcn, iAC), 1 : nAC);           
                obj.LogFcn('Vehicle (aircraft) data extracted!');
                
                %Extract data for the profile(s)
                obj.LogFcn('Extracting the profile information...');
                %Fuselages ...
                ProfileData.Fuselage = obj.extractProfileData(obj.FuselageProfiles);
                obj.LogFcn(sprintf('%i fuselage profiles extracted...', numel(ProfileData.Fuselage)));
                %Aerofoils ...
                ProfileData.Wing     = obj.extractProfileData(obj.WingProfiles);
                obj.LogFcn(sprintf('%i aerofoil profiles extracted...', numel(ProfileData.Wing)));
                obj.LogFcn('Profile library generated!');
                
                %Consolidate the aircraft & profile data
                %   - Link profiles to the different beam elements
                obj.LogFcn('Matching profiles to wing/fuselage elements...');
                Aircraft = arrayfun(@(iAC) obj.consolidateVehicleData( ...
                    AircraftData(iAC), ProfileData, obj.LogFcn, iAC), 1 :nAC);
                obj.LogFcn('%i Aircraft successfully consolidated!');
                
                %If we get this far then we have finished the data
                %extraction without error. Exit while loop!
                obj.EndConversion = true;
                
            end
            
            %Did we invoke any errors during the conversion?
            if obj.QuitConversion
                return
            end               

            
            %% Generate the AWI Framework model using the CPACS data
            AWIObj = obj.generateAWIobjects(Aircraft, obj.LogFcn);
            
        end
        function exportAWIToCPACs(obj, AWIFramework)
            
            
        end        
    end
    
    methods (Access = private) %Error handling
        function terminateConversion(obj, token)
            %terminateConversion Halts the CPACS-to-AWI conversion process
            %and informs the user of the reason why.
            
            if nargin < 2
                token = '-UNKNOWN TOKEN-';
            end
            
            %Format message
            str = sprintf(['Could not find the token ''%s'' in the ', ...
                'CPACS .xml file, unable to proceed. Quitting ', ...
                'CPACS conversion.'], token);
            
            %Tell the user
            obj.LogFcn(str);
            
            %Tell the object
            obj.EndConversion  = true;
            obj.QuitConversion = true;
            
        end
    end
    
    methods (Static, Access = private) %Extracting data from the master CPACS object
        
        %extractProfileData
        function ProfileData = extractProfileData(pf_Data)
            %extractProfileData Extracts the profile data for the fuselage
            %and aerofoils in the collection.
            %
            % 'pf_Data' is a collector for all the profiles in the
            % model. Each profile is a child of 'pf_Data' and is an
            % 'awi.model.Entity' due to the simplistic nature of the
            % AWI .xml importer.
            %
            % This code has been tested for CPACS tokens
            % "fuselageProfile" & "wingAirfoil" and both work without
            % error.
            
            %Check for empty
            %   - Can get an empty 'awi.model.Entity' if 'findall' does
            %   not find anything in the collection.
            if isempty(pf_Data)
                ProfileData = [];
                return
            end
            
            %Define the fields of the structure array
            ProfileData = struct( ...
                'uID'        , [], ...
                'Name'       , [], ...
                'Description', [], ...
                'PointList'  , []);
            
            %ASSUME all children of the collector are valid fuselage
            %profiles
            nFuseProfile = numel(pf_Data.Children);
            
            %Make sure we have the correct amount of structures!
            ProfileData = repmat(ProfileData, [1, nFuseProfile]);
            
            %Assign the "Name" and "Description" fields
            nam = get(pf_Data.Children, 'Name');
            des = get(pf_Data.Children, 'Description');
            des = vertcat(des{:}); %ASSUME that des returns a cell array of cell arrays...
            [ProfileData(:).Name] = deal(nam{:});
            [ProfileData(:).Description] = deal(des{:});
            
            %Grab the points
            for iP = 1 : nFuseProfile
                %Find the "pointList"
                pl = findall(pf_Data.Children(1), 'Name', 'pointList');
                if isempty(pl)
                    continue
                end
                %Extract the data
                ProfileData(iP).PointList = [ ...
                    str2num(pl.x), ...
                    str2num(pl.y), ...
                    str2num(pl.z)]; %#ok<ST2NM>
            end
            
        end
        
        %extractAircraftData
        function AircraftData = extractAircraftData(ac_Data, logfcn, acx)
            %extractAircraftData 
            
            if nargin < 2
               logfcn = @(s) fprintf('%s\n', s); 
            end
            if nargin < 3
                acx = 1;
            end
            
            %We are currently one level below the master CPACS element
            lvl = 1;
            
            logfcn(sprintf('Extracting data for aircraft %i', acx));
            
            %Start with something sensible
            AircraftData = struct();
            
            %Assume that the first child is the "model" element of the
            %aircraft
            mdl = ac_Data.Children(1);            
            if isempty(mdl)
                return
            end
            
            %Get reference data
            logfcn([blanks(lvl), 'Extracting reference data'])
            AircraftData.Reference = i_extractReferenceData(mdl);
            
            %Get handles to all fuselages
            f = findall(ac_Data, 'Name', 'fuselages');
            if isempty(f)
                logfcn([blanks(lvl), '* * No fuselage data in the model * *']);
                AircraftData.Fuselages = [];
            else
                logfcn([blanks(lvl), 'Extracting fuselage data']);
                %ASSUME all children are valid beam/section objects...
                AircraftData.Fuselages = arrayfun(@(iF) i_extractBeamElementData( ...
                    f.Children(iF), logfcn, 'fuselage', iF), 1 : numel(f.Children));
            end
            
            %Get handles to all wings
            w = findall(ac_Data, 'Name', 'wings');
            if isempty(w)
                logfcn([blanks(lvl), '* * No wing data in the model * *']);
                AircraftData.Wings = [];
            else
                logfcn([blanks(lvl), 'Extracting wing data']);
                %ASSUME all children are valid beam/section objects...
                AircraftData.Wings = arrayfun(@(iW) i_extractBeamElementData( ...
                    w.Children(iW), logfcn, 'wing', iW), 1 : numel(w.Children));
            end
            
            %Extract reference data
            function Reference = i_extractReferenceData(mdl)
                %i_extractReferenceData Extracts the reference data from
                %the CPACS "model" element. 
                %
                % Searches for the following data:
                %   - 'length'
                %   - 'area'
                %   - 'point'
                
                %Start with an empty structure
                Reference = struct();
                
                %Find the "reference" element
                ref = findall(mdl, 'Name', 'reference');                
                if isempty(ref)
                    Reference = [];
                    return
                end
                
                %Length?
                if isprop(ref, 'length')
                    Reference.Length = str2double(ref.length);
                end
                
                %Area?
                if isprop(ref, 'area')
                    Reference.Area = str2double(ref.area);
                end
                
                %Point?
                point = findall(ref, 'Name', 'point');
                if ~isempty(point)
                    Reference.Point = [ ...
                        str2double(point.x) ; ...
                        str2double(point.y) ; ...
                        str2double(point.z)];
                end
                
                %Check if we have actually populated the structure
                if numel(fieldnames(Reference)) == 0
                    Reference = [];
                end
                
            end
            
            %Extract beam data
            function BeamData = i_extractBeamElementData(b_data, logfcn, type, bmx)
                %i_extractBeamElementData Extracts data from a CPACS beam
                %element.
                %
                % A CPACS beam element is assumed to comprise of:
                %   - Transformation(s)
                %   - Section(s)
                %   - Positioning(s)
                %   - Segment(s)
                
                indent = 3;
                
                %Template for a blank CPACS "transformation" element
                T = struct( ...
                    'Scaling'    , [0, 0, 0], ...
                    'Rotation'   , [0, 0, 0] , ...
                    'Translation', [0, 0, 0]);
                
                %Start with the basic format of a CPACS beam
                BeamData = struct( ... 
                    'Transformation', T, ...
                    'Sections'      , [], ...
                    'Positionings'  , [], ...
                    'Segments'      , []);
                
                logfcn([blanks(indent-1), sprintf('Extracting data for %s element %i', type, bmx)]);
                
                if isprop(b_data, 'name')        %Grab name        
                    BeamData.Name = b_data.name;
                else
                    BeamData.Name = '';
                end
                if isprop(b_data, 'description') %Grab description 
                    BeamData.Description = b_data.description;
                else
                    BeamData.Description = '';
                end
                
                %Extract beam transformation data
                t = find(b_data.Children, 'Name', 'transformation');
                if isempty(t)
                    logfcn([blanks(indent), '* * Default transformation data used * *']);
                else
                    logfcn([blanks(indent), 'Extracting transformation data.']);
                    BeamData.Transformation = i_extractTransform(t);
                    logfcn([blanks(indent), 'Transformation data extracted!']);
                end
                   
                %Extract sections
                s = find(b_data.Children, 'Name', 'sections');
                if isempty(s)
                    logfcn([blanks(indent), '* * No section data found * *']);
                else
                    logfcn([blanks(indent), 'Extracting section data...']);
                    BeamData.Sections = i_extractSections(s);
                    logfcn([blanks(indent), sprintf('%i sections extracted!', ...
                        numel(BeamData.Sections))]);
                end
                
                %Extract positionings
                p = find(b_data.Children, 'Name', 'positionings');
                if isempty(p)
                    logfcn([blanks(indent), '* * No position data found * *']);
                else
                    logfcn([blanks(indent), 'Extracting position data...']);
                    BeamData.Positionings = i_extractPositionings(p);
                    logfcn([blanks(indent), sprintf('%i position data sets extracted!', ...
                       numel(BeamData.Positionings))]);
                end
                
                %Extract segments
                sg = find(b_data.Children, 'Name', 'segments');
                if isempty(p)
                    logfcn([blanks(indent), '* * No segment data found * *']);
                else
                    logfcn([blanks(indent), 'Extracting segment data...']);
                    BeamData.Segments = i_extractSegments(sg);
                    logfcn([blanks(indent), sprintf('%i segment data sets extracted!', ...
                       numel(BeamData.Segments))]);
                end     
                
                logfcn([blanks(indent-1), sprintf('Data extracted for %s element %i!', type, bmx)]);
                
                
                %...Extract "transformation" element
                function Transform = i_extractTransform(trnsfrm)
                    %i_extractTransform Extracts transformation data.
                    %
                    %   - A transformation has 3 fields: Translation,
                    %   Rotation & Scaling.
                    %   - The defualt value for each field is [0,0,0] 
                    
                    %Translation
                    t_ = find(trnsfrm.Children, 'Name', 'translation');
                    if isempty(t_) %Use default or extract data ...
                        Transform.Translation = T.Translation;
                    else
                        Transform.Translation = [ ...
                            str2double(t_.x), ...
                            str2double(t_.y), ...
                            str2double(t_.z)];
                    end
                    
                    %Rotation
                    r_ = find(trnsfrm.Children, 'Name', 'rotation');
                    if isempty(r_) %Use default or extract data ...
                        Transform.Rotation = T.Rotation;
                    else
                        Transform.Rotation = [ ...
                            str2double(r_.x), ...
                            str2double(r_.y), ...
                            str2double(r_.z)];
                    end
                    
                    %Scaling
                    s_ = find(trnsfrm.Children, 'Name', 'scaling');
                    if isempty(s_) %Use default or extract data ...
                        Transform.Scaling = T.Rotation;
                    else
                        Transform.Scaling = [ ...
                            str2double(s_.x), ...
                            str2double(s_.y), ...
                            str2double(s_.z)];
                    end
                    
                end
                
                %...Extract "section" element(s)
                function Sections = i_extractSections(sects)
                    %i_extractSections Extracts the section data.
                    %
                    %   - A section has the following fields: Name,
                    %   Description, Transformation, Elements
                    
                    %Start with default structure
                    Sections = struct( ...
                        'Name'          , [], ...
                        'Description'   , [], ...
                        'Transformation', [], ...
                        'Elements'      , []);
                    Elements = struct( ...
                        'Name'          , [], ...
                        'Description'   , [], ...
                        'ProfileUID'    , [], ...
                        'Transformation', []);
                    
                    %ASSUME all children of 'sects' are valid section
                    %elements.
                    nSect = numel(sects.Children);
                    
                    %Preallocate
                    Sections = repmat(Sections, [1, nSect]);
                    
                    nam = get(sects.Children, 'Name');
                    des = get(sects.Children, 'Description');
                    [Sections(:).Name]        = deal(nam{:});
                    [Sections(:).Description] = deal(des{:});
                    
                    %Extra transformation and elements data individually
                    for iS = 1 : nSect
                        
                        %...Transformation
                        t_ = find(sects.Children(iS), 'Name', 'transformation');
                        if isempty(t_)
                            Sections(iS).Transformation = T;
                        else
                            Sections(iS).Transformation = i_extractTransform(t_);
                        end     
                        
                        %...Elements
                        elms = find(sects.Children(iS).Children, 'Name', 'elements');
                        if isempty(elms) %Go to next section if no elements are found
                            continue
                        end
                        
                        %Go one level down from the collector
                        elms  = elms.Children;
                        nElem = numel(elms);
                        Elem  = repmat(Elements, [1, nElem]);
                        
                        %Grab and assign the name, description & profile ID
                        nam    = {elms.Name};
                        des    = {elms.Description};
                        %prfUID = {elms.profileUID};
                        [Elem(:).Name]        = deal(nam{:});
                        [Elem(:).Description] = deal(des{:});
%                         [Elem(:).ProfileUID]  = deal(prfUID{:});
                        
                        %Grab "profileUID" (if applicable)
                        idx = arrayfun(@(p) isprop(p, 'profileUID'), elms);
                        val = get(elms(idx), 'profileUID');
                        if ~iscell(val) %Force cell
                            val = {val}; 
                        end
                        [Elem(idx).ProfileUID] = deal(val{:});                    
                        
                        %Grab "airfoilUID" (if applicable)
                        idx = arrayfun(@(p) isprop(p, 'airfoilUID'), elms);
                        val = get(elms(idx), 'airfoilUID');
                        if ~iscell(val) %Force cell
                            val = {val}; 
                        end
                        [Elem(idx).ProfileUID] = deal(val{:}); 
                        
                        %Extract each element's transformation data
                        for iE = 1 : nElem
                            t_ = find(elms.Children, 'Name', 'transformation');
                            if isempty(t_)
                                Elem(iE).Transformation = T;
                            else
                                Elem(iE).Transformation = i_extractTransform(t_);
                            end
                        end
                        
                        %Assign Elements data to the Sections structure
                        Sections(iS).Elements = Elem;
                        
                    end
                    
                end
                
                %...Extract "positioning" element(s)
                function Positionings = i_extractPositionings(pos)
                    %i_extractPositionings Extracts the position data.
                    %
                    %   - A position has the following fields: Name,
                    %   Length, SweepAngle, DihedralAngle, FromSectionUID, 
                    %   ToSectionUID
                    
                    %Start with default structure
                    Positionings = struct( ...
                        'Name'          , [], ...
                        'Length'        , [], ...
                        'SweepAngle'    , [], ...
                        'DihedralAngle' , [], ...
                        'FromSectionUID', [], ...
                        'ToSectionUID'  , []);
                    
                    %ASSUME all children of 'pos' are valid positioning
                    %elements.
                    nPos = numel(pos.Children);
                    
                    %Preallocate
                    Positionings = repmat(Positionings, [1, nPos]);
                    
                    %Grab "names"
                    idx = arrayfun(@(p) isprop(p, 'name'), pos.Children);
                    val = get(pos.Children(idx), 'name');
                    [Positionings(idx).Name] = deal(val{:});
                    
                    %Grab "length"
                    idx = arrayfun(@(p) isprop(p, 'length'), pos.Children);
                    val = get(pos.Children(idx), 'length');
                    val = num2cell(str2double(val)); %Assume length is a scalar
                    [Positionings(idx).Length] = deal(val{:});
                    
                    %Grab "sweepAngle"
                    idx = arrayfun(@(p) isprop(p, 'sweepAngle'), pos.Children);
                    val = get(pos.Children(idx), 'sweepAngle');
                    val = num2cell(str2double(val)); %Assume sweep angle is a scalar
                    [Positionings(idx).SweepAngle] = deal(val{:});
                                        
                    %Grab "dihedralAngle"
                    idx = arrayfun(@(p) isprop(p, 'dihedralAngle'), pos.Children);
                    val = get(pos.Children(idx), 'dihedralAngle');
                    val = num2cell(str2double(val)); %Assume dihedral angle is a scalar
                    [Positionings(idx).DihedralAngle] = deal(val{:});
                    
                    %Grab "fromSectionUID"
                    idx = arrayfun(@(p) isprop(p, 'fromSectionUID'), pos.Children);
                    val = get(pos.Children(idx), 'fromSectionUID');
                    if ~iscell(val) %Force cell
                        val = {val};
                    end
                    [Positionings(idx).FromSectionUID] = deal(val{:});
                    
                    %Grab "toSectionUID"
                    idx = arrayfun(@(p) isprop(p, 'toSectionUID'), pos.Children);
                    val = get(pos.Children(idx), 'toSectionUID');
                    if ~iscell(val) %Force cell
                        val = {val};
                    end
                    [Positionings(idx).ToSectionUID] = deal(val{:});
                    
                end
                
                %...Extract "segment" element(s)
                function Segments = i_extractSegments(seg)
                    %i_extractSegments Extracts the segment data.
                    %
                    %   - A position has the following fields: Name,
                    %   Description, FromElementUID, ToElementUID
                    
                    %Start with default structure
                    Segments = struct( ...
                        'Name'          , [], ...
                        'Description'   , [], ...
                        'FromElementUID', [], ...
                        'ToElementUID'  , []);
                    
                    %ASSUME all children of 'pos' are valid positioning
                    %elements.
                    nSeg = numel(seg.Children);
                    
                    %Preallocate
                    Segments = repmat(Segments, [1, nSeg]);
                    
                    %Grab "names"
                    idx = arrayfun(@(p) isprop(p, 'name'), seg.Children);
                    val = get(seg.Children(idx), 'name');
                    if ~iscell(val) %Force cell
                        val = {val};
                    end
                    [Segments(idx).Name] = deal(val{:});
                    
                    %Grab "description"
                    idx = arrayfun(@(p) isprop(p, 'description'), seg.Children);
                    val = get(seg.Children(idx), 'description');
                    if ~iscell(val) %Force cell
                        val = {val};
                    end
                    [Segments(idx).Description] = deal(val{:});
                    
                    %Grab "fromSectionUID"
                    idx = arrayfun(@(p) isprop(p, 'fromElementUID'), seg.Children);
                    val = get(seg.Children(idx), 'fromElementUID');
                    if ~iscell(val) %Force cell
                        val = {val};
                    end
                    [Segments(idx).FromElementUID] = deal(val{:});
                    
                    %Grab "toSectionUID"
                    idx = arrayfun(@(p) isprop(p, 'toElementUID'), seg.Children);
                    val = get(seg.Children(idx), 'toElementUID');
                    if ~iscell(val) %Force cell
                        val = {val};
                    end
                    [Segments(idx).ToElementUID] = deal(val{:});
                    
                end
                
            end
                        
        end

        %consolidateAircraftData 
        function Aircraft = consolidateVehicleData(Aircraft, Profiles, logfcn, acx)
            
            %Wing(s)
            
            %Get the profile names for indexing
            pf_names = {Profiles.Wing.Name};
            
            
            
            %Fuselage(s)
            
            
            
        end
        
        %generateAWIobjects
        function AWI_Models = generateAWIobjects(Aircraft, logfcn)
            
            
            
        end
        
    end
    
    
    
    
end

