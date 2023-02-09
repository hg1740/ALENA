function output = sortFameInput(FAME,fcn)

if nargin < 2
    fcn = @(s) fprintf('%s\n',s);
end

% Extract the fieldnames 
User_var = fieldnames(FAME.USER_VARIABLES);

% Structure in which the user variables will be stored
output.UserVariables = [];

for i = 1:length(User_var)

    % Remove the 'eval' and '&' strings from the input line
    var = eval(regexprep(FAME.USER_VARIABLES.(User_var{i,1}){1,1},{'&','eval'},''));
    
    % Use the eval function to make the string a variable and assign var
    eval([User_var{i,1} ' = var;']);
    
    % Store these values in a structure for future reference
    output.UserVariables.(User_var{i,1}) = eval([User_var{i,1} ';']);
end

%% Wing Reference point
ref_point = [eval(regexprep(FAME.AIRPLANE.GLOBAL_DATA.WING_POSITION.X_REF{1},{'&','eval'},'')),...
    0,eval(regexprep(FAME.AIRPLANE.GLOBAL_DATA.WING_POSITION.Z_REF{1},{'&','eval'},''))]/1000;

%% Dihedral
dihedral = zeros(1,2);
for i = 1:length(FAME.WING.GLOBAL_GEOMETRY.DIHEDRAL_DISTRIBUTION.DIHEDRAL)
    dihedral(i,1) = eval(regexprep(FAME.WING.GLOBAL_GEOMETRY.DIHEDRAL_DISTRIBUTION.ETA_BEG{i},{'&','eval'},''));
    dihedral(i,2) = eval(regexprep(FAME.WING.GLOBAL_GEOMETRY.DIHEDRAL_DISTRIBUTION.DIHEDRAL{i},{'&','eval'},''));
end
dihedral(end+1,:) = [1,dihedral(end,2)];

%% Chord Distribution
chord = zeros(1,2);
for i = 1:length(FAME.WING.GLOBAL_GEOMETRY.CHORD_DISTRIBUTION.CHORD)
    chord(i,1) = eval(regexprep(FAME.WING.GLOBAL_GEOMETRY.CHORD_DISTRIBUTION.ETA{i},{'&','eval'},''));
    chord(i,2) = eval(regexprep(FAME.WING.GLOBAL_GEOMETRY.CHORD_DISTRIBUTION.CHORD{i},{'&','eval'},''))/1000;
end

%% Quarter Chord Sweep
if isfield(FAME.WING.GLOBAL_GEOMETRY.SWEEP_DISTRIBUTION,'SWEEP_LE')
    SweepType = 'SWEEP_LE';
    SweepLoc  = 0.0;
elseif isfield(FAME.WING.GLOBAL_GEOMETRY.SWEEP_DISTRIBUTION,'SWEEP_25')
    SweepType = 'SWEEP_25';
    SweepLoc  = 0.25;
elseif isfield(FAME.WING.GLOBAL_GEOMETRY.SWEEP_DISTRIBUTION,'SWEEP_TE')
    SweepType = 'SWEEP_TE';
    SweepLoc  = 1.0;
else
    error('Neither SWEEP_LE nor SWEEP_25 have been found for the wing');
end

QCsweep = zeros(1,2);
for i = 1:length(FAME.WING.GLOBAL_GEOMETRY.SWEEP_DISTRIBUTION.(SweepType))
    QCsweep(i,1) = eval(regexprep(FAME.WING.GLOBAL_GEOMETRY.SWEEP_DISTRIBUTION.ETA_BEG{i},{'&','eval'},''));
    QCsweep(i,2) = eval(regexprep(FAME.WING.GLOBAL_GEOMETRY.SWEEP_DISTRIBUTION.(SweepType){i},{'&','eval'},''));
end
QCsweep(end+1,:) = [1,QCsweep(end,2)];
QCsweep(:,3)     = SweepLoc;

%% Front Spar Position
fspar = zeros(1,2);
for i = 1:length(FAME.WING.DETAILED_GEOMETRY.FRONT_SPAR_POS.ETA)
    fspar(i,1) = eval(regexprep(FAME.WING.DETAILED_GEOMETRY.FRONT_SPAR_POS.ETA{i},{'&','eval'},''));
    fspar(i,2) = eval(regexprep(FAME.WING.DETAILED_GEOMETRY.FRONT_SPAR_POS.CHORD_REL{i},{'&','eval'},''));
end

%% Rear Spar Position
aspar = zeros(1,2);
for i = 1:length(FAME.WING.DETAILED_GEOMETRY.REAR_SPAR_POS.ETA)
    aspar(i,1) = eval(regexprep(FAME.WING.DETAILED_GEOMETRY.REAR_SPAR_POS.ETA{i},{'&','eval'},''));
    aspar(i,2) = eval(regexprep(FAME.WING.DETAILED_GEOMETRY.REAR_SPAR_POS.CHORD_REL{i},{'&','eval'},''));
end

%% Thickness to Chord Ratio
thickness = zeros(1,2);
for i = 1:length(FAME.WING.GLOBAL_GEOMETRY.THICKNESS_DISTRIBUTION.THICKNESS_REL)
    thickness(i,1) = eval(regexprep(FAME.WING.GLOBAL_GEOMETRY.THICKNESS_DISTRIBUTION.ETA{i},{'&','eval'},''));
    thickness(i,2) = eval(regexprep(FAME.WING.GLOBAL_GEOMETRY.THICKNESS_DISTRIBUTION.THICKNESS_REL{i},{'&','eval'},''));
end

if ~isfield(FAME.USER_VARIABLES,'rib_spc')
    fprintf('\nWARNING - NO RIB SPACING FOUND. ASSIGNED VALUE OF RP = 0.65m');
    rib_pitch = 0.65; %
else
    rib_pitch = eval(regexprep(FAME.USER_VARIABLES.rib_spc{1,1},{'&','eval'},''))/1000;
end

if ~isfield(FAME.AIRPLANE,'HORIZONTAL_TAILPLANE_POSITION')
    
    if isfield(FAME.HORIZONTAL_TAILPLANE.GLOBAL_GEOMETRY.SWEEP_DISTRIBUTION,'SWEEP_LE')
        SweepType = 'SWEEP_LE';
        SweepLoc  = 0.0;
    elseif isfield(FAME.HORIZONTAL_TAILPLANE.GLOBAL_GEOMETRY.SWEEP_DISTRIBUTION,'SWEEP_25')
        SweepType = 'SWEEP_25';
        SweepLoc  = 0.25;
    else
        error('Neither SWEEP_LE nor SWEEP_25 have been found for the wing');
    end
    % Quarter Chord Sweep
    HTPQCsweep = zeros(1,2);
    for i = 1:length(FAME.HORIZONTAL_TAILPLANE.GLOBAL_GEOMETRY.SWEEP_DISTRIBUTION.(SweepType))
        HTPQCsweep(i,1) = eval(regexprep(FAME.HORIZONTAL_TAILPLANE.GLOBAL_GEOMETRY.SWEEP_DISTRIBUTION.ETA_BEG{i},{'&','eval'},''));
        HTPQCsweep(i,2) = eval(regexprep(FAME.HORIZONTAL_TAILPLANE.GLOBAL_GEOMETRY.SWEEP_DISTRIBUTION.(SweepType){i},{'&','eval'},''));
    end
    HTPQCsweep(end+1,:) = [1,HTPQCsweep(end,2)];
    HTPQCsweep(:,3)     = SweepLoc;
    
    % Thickness to Chord Ratio
    HTPthickness = zeros(1,2);
    for i = 1:length(FAME.HORIZONTAL_TAILPLANE.GLOBAL_GEOMETRY.THICKNESS_DISTRIBUTION.THICKNESS_REL)
        HTPthickness(i,1) = eval(regexprep(FAME.HORIZONTAL_TAILPLANE.GLOBAL_GEOMETRY.THICKNESS_DISTRIBUTION.ETA{i},{'&','eval'},''));
        HTPthickness(i,2) = eval(regexprep(FAME.HORIZONTAL_TAILPLANE.GLOBAL_GEOMETRY.THICKNESS_DISTRIBUTION.THICKNESS_REL{i},{'&','eval'},''));
    end
    
    % Dihedral
    HTPdihedral = zeros(1,2);
    for i = 1:length(FAME.HORIZONTAL_TAILPLANE.GLOBAL_GEOMETRY.DIHEDRAL_DISTRIBUTION.DIHEDRAL)
        HTPdihedral(i,1) = eval(regexprep(FAME.HORIZONTAL_TAILPLANE.GLOBAL_GEOMETRY.DIHEDRAL_DISTRIBUTION.ETA_BEG{i},{'&','eval'},''));
        HTPdihedral(i,2) = eval(regexprep(FAME.HORIZONTAL_TAILPLANE.GLOBAL_GEOMETRY.DIHEDRAL_DISTRIBUTION.DIHEDRAL{i},{'&','eval'},''));
    end
    HTPdihedral(end+1,:) = [1,HTPdihedral(end,2)];
    
    % Chord Distribution
    HTPchord =zeros(1,2);
    for i = 1:length(FAME.HORIZONTAL_TAILPLANE.GLOBAL_GEOMETRY.CHORD_DISTRIBUTION.CHORD)
        HTPchord(i,1) = eval(regexprep(FAME.HORIZONTAL_TAILPLANE.GLOBAL_GEOMETRY.CHORD_DISTRIBUTION.ETA{i},{'&','eval'},''));
        HTPchord(i,2) = eval(regexprep(FAME.HORIZONTAL_TAILPLANE.GLOBAL_GEOMETRY.CHORD_DISTRIBUTION.CHORD{i},{'&','eval'},''))/1000;
    end
    
    HTPspan = eval(regexprep(FAME.HORIZONTAL_TAILPLANE.GLOBAL_GEOMETRY.SPAN{1},{'&','eval'},''))/1000;

    if isfield(FAME.AIRPLANE.GLOBAL_DATA,'HORIZONTAL_TAILPLANE_POSITION')
        if isfield(FAME.AIRPLANE.GLOBAL_DATA.HORIZONTAL_TAILPLANE_POSITION,'LEVER_ARM_X')
             HTPref_point = [ref_point(1),0,0] + [eval(regexprep(FAME.AIRPLANE.GLOBAL_DATA.HORIZONTAL_TAILPLANE_POSITION.LEVER_ARM_X{1},{'&','eval'},'')),...
                 0,eval(regexprep(FAME.AIRPLANE.GLOBAL_DATA.HORIZONTAL_TAILPLANE_POSITION.Z_REF{1},{'&','eval'},''))]/1000;
        else
            HTPref_point = [eval(regexprep(FAME.AIRPLANE.GLOBAL_DATA.HORIZONTAL_TAILPLANE_POSITION.X_REF{1},{'&','eval'},'')),...
                 0,eval(regexprep(FAME.AIRPLANE.GLOBAL_DATA.HORIZONTAL_TAILPLANE_POSITION.Z_REF{1},{'&','eval'},''))]/1000;
        end
    else
        fprintf('\nWARNING - NO HTP POSITION DATA!!!');
        HTPref_point = [];
    end
    if ~isfield(FAME.USER_VARIABLES,'rib_spc')
        fprintf('\nWARNING - NO RIB SPACING FOUND. ASSIGNED VALUE OF RP = 0.65m');
        HTPrib_pitch = 0.65; %
    else
        HTPrib_pitch = eval(regexprep(FAME.USER_VARIABLES.rib_spc{1,1},{'&','eval'},''))/1000;
    end
    
else
    HTPQCsweep      = [];
    HTPchord        = [];
    HTPdihedral     = [];
    HTPthickness    = [];
    HTPspan         = [];
    HTPref_point    = [];
    HTPrib_pitch    = [];
end


%% Weights
MTOW = eval(regexprep(FAME.AIRPLANE.GLOBAL_DATA.DESIGN_WEIGHTS.MTOW{1},{'&','eval'},''));
MZFW = eval(regexprep(FAME.AIRPLANE.GLOBAL_DATA.DESIGN_WEIGHTS.MZFW{1},{'&','eval'},''));
MLW  = eval(regexprep(FAME.AIRPLANE.GLOBAL_DATA.DESIGN_WEIGHTS.MLW{1},{'&','eval'},''));

%% Engines
Engine = [];
if isfield(FAME,'PROPULSION')
    for i = 1:numel(FAME.PROPULSION.ID)
        Engine(i).eta  = eval(regexprep(FAME.PROPULSION.ETA{1},{'&','eval'},''));
        Engine(i).engine_mass = eval(regexprep(FAME.PROPULSION.ENGINE_MASS{1},{'&','eval'},''));
        if isfield(FAME.PROPULSION,'AXIS_Z_LOC')
            Engine(i).engine_z_loc = eval(regexprep(FAME.PROPULSION.AXIS_Z_LOC{1},{'&','eval'},''))/1000;
        else
            Engine(i).engine_z_loc = eval(regexprep(FAME.PROPULSION.POD_PYLON_CG_Z_LOC{1},{'&','eval'},''))/1000;
        end
        Engine(i).pylon_mass = eval(regexprep(FAME.PROPULSION.PYLON_MASS{1},{'&','eval'},''));
        Engine(i).pylon_x_loc = eval(regexprep(FAME.PROPULSION.POD_PYLON_CG_X_LOC{1},{'&','eval'},''))/1000;
        Engine(i).pylon_z_loc = eval(regexprep(FAME.PROPULSION.POD_PYLON_CG_Z_LOC{1},{'&','eval'},''))/1000;
        
    end
end

%% Control Surface Information

% Perform a search on the control surfaces: we assume that the control
% surface has the identifier geo in its name

% Get all the fieldnames
geonames = fieldnames(FAME.WING.DETAILED_GEOMETRY);

% Find the control surfaces
csnameidx = ~cellfun(@isempty,regexpi(geonames,'geo'));

% Get the control surface headers
csnames = geonames(csnameidx);

for i = 1:numel(csnames)
    
    % Remove _GEO from the header 
    csname    = lower(regexprep(csnames{i},{'_','GEO'},''));
    
    % Get the control surface properties
    csprops   = fieldnames(FAME.WING.DETAILED_GEOMETRY.(csnames{i}));

    % Assign the properties to the structure
    for j = 1:numel(csprops)
        Wing.ControlSurfaces.(csname).(lower(csprops{j})) = [];
    end
    
    % Iterate through all the eta values of that control surface
    for k = 1:length(FAME.WING.DETAILED_GEOMETRY.(csnames{i}).(csprops{j}))
        
        % Get the properties of the control surface
        for j = 1:numel(csprops)
            prop  = eval(regexprep(FAME.WING.DETAILED_GEOMETRY.(csnames{i}).(csprops{j}){k},{'&','eval'},''));
            Wing.ControlSurfaces.(csname).(lower(csprops{j}))  = [Wing.ControlSurfaces.(csname).(lower(csprops{j}));prop];
        end
    end
    
end

Wing.dihedral       = dihedral;
Wing.sweep          = QCsweep;
Wing.thickness      = thickness ;
Wing.chord          = chord;
Wing.span           = eval(regexprep(FAME.WING.GLOBAL_GEOMETRY.SPAN{1},{'&','eval'},''))/1000;
Wing.ref_point      = ref_point;
Wing.fspar          = fspar;
Wing.aspar          = aspar;
Wing.rib_pitch      = rib_pitch;
Wing.stringer_pitch = 0.2; % TODO - Hardcoded stringer pitch!!
Wing.Material.Setup      = [];
Wing.Material.Allowables = [];

HTP.dihedral       = HTPdihedral;
HTP.sweep          = HTPQCsweep;
HTP.thickness      = HTPthickness ;
HTP.chord          = HTPchord;
HTP.span           = HTPspan;
HTP.ref_point      = HTPref_point;
HTP.rib_pitch      = HTPrib_pitch;

Weights.mtow    = MTOW;
Weights.mzfw    = MZFW;
Weights.mlw     = MLW;

output.Wing     = Wing;
output.HTP      = HTP;
output.Weights  = Weights;
output.Engine   = Engine;


end

