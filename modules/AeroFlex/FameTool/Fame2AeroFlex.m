% The following function is in charge of converting a FAME structure,
% derived directly from the *.fm4 file, into one that is useable by
% AeroFlex. The purpose of this script is not to process the data as such,
% but more to import the data. The data is then interpolated by the setup
% box function.

function param = Fame2AeroFlex(FAME)

%% WING PLANFORM PROPERTIES

%% EVALUATE USER VARIABLES

User_var = fieldnames(FAME.USER_VARIABLES);

for i = 1:length(User_var)
    
    % TODO Evaluate whether there is an evaluate function in the variable
    eval([User_var{i,1} '= str2double(FAME.USER_VARIABLES.' User_var{i,1} '{1,1});']);
end

%% Wing Reference point
ref_point = [str2double(FAME.AIRPLANE.GLOBAL_DATA.WING_POSITION.X_REF{1})/1000,...
    0,str2double(FAME.AIRPLANE.GLOBAL_DATA.WING_POSITION.Z_REF{1})/1000];

%% Extract some key wing properties from the FAME file

% Quarter Chord Sweep
QCsweep = [];
for i = 1:length(FAME.WING.GLOBAL_GEOMETRY.SWEEP_DISTRIBUTION.SWEEP_25)
    
    check_var = strtrim(regexp(FAME.WING.GLOBAL_GEOMETRY.SWEEP_DISTRIBUTION.ETA_BEG{i},'(?<=&).*()','match'));
    
    Eta = str2double(FAME.WING.GLOBAL_GEOMETRY.SWEEP_DISTRIBUTION.ETA_BEG{i});
    % If the input is defined by a user variable then call it
    if ~isempty(check_var)
        eval(['Eta = ' check_var{1} ';']);
    end
    
    QCsweep = [QCsweep; [Eta,...
        str2double(FAME.WING.GLOBAL_GEOMETRY.SWEEP_DISTRIBUTION.SWEEP_25{i})]];
    
end
QCsweep(end+1,:) = [1,QCsweep(end,2)];

% Thickness to Chord Ratio
thickness = [];
for i = 1:length(FAME.WING.GLOBAL_GEOMETRY.THICKNESS_DISTRIBUTION.THICKNESS_REL)
    
    check_var = strtrim(regexp(FAME.WING.GLOBAL_GEOMETRY.THICKNESS_DISTRIBUTION.ETA{i},'(?<=&).*()','match'));
    
    Eta = str2double(FAME.WING.GLOBAL_GEOMETRY.THICKNESS_DISTRIBUTION.ETA{i});
    
    % If the input is defined by a user variable then call it
    if ~isempty(check_var)
        eval(['Eta = ' check_var{1} ';']);
    end
    
    thickness = [thickness; [Eta,...
        str2double(FAME.WING.GLOBAL_GEOMETRY.THICKNESS_DISTRIBUTION.THICKNESS_REL{i})]];
    
end

% Dihedral
dihedral = [];
for i = 1:length(FAME.WING.GLOBAL_GEOMETRY.DIHEDRAL_DISTRIBUTION.DIHEDRAL)
    
    check_var = strtrim(regexp(FAME.WING.GLOBAL_GEOMETRY.DIHEDRAL_DISTRIBUTION.ETA_BEG{i},'(?<=&).*()','match'));
    
    Eta = str2double(FAME.WING.GLOBAL_GEOMETRY.CHORD_DISTRIBUTION.ETA{i});
    
    % If the input is defined by a user variable then call it
    if ~isempty(check_var)
        eval(['Eta = ' check_var{1} ';']);
    end
    
    dihedral = [dihedral; [Eta,...
        str2double(FAME.WING.GLOBAL_GEOMETRY.DIHEDRAL_DISTRIBUTION.DIHEDRAL{i})]];
    
end
dihedral(end+1,:) = [1,dihedral(end,2)];

% Chord Distribution
chord = [];
for i = 1:length(FAME.WING.GLOBAL_GEOMETRY.CHORD_DISTRIBUTION.CHORD)
    
    check_var = strtrim(regexp(FAME.WING.GLOBAL_GEOMETRY.CHORD_DISTRIBUTION.ETA{i},'(?<=&).*()','match'));
    
    Eta = str2double(FAME.WING.GLOBAL_GEOMETRY.CHORD_DISTRIBUTION.ETA{i});
    
    % If the input is defined by a user variable then call it
    if ~isempty(check_var)
        eval(['Eta = ' check_var{1} ';']);
    end
    
    chord = [chord; [Eta,...
        str2double(FAME.WING.GLOBAL_GEOMETRY.CHORD_DISTRIBUTION.CHORD{i})/1000]];
    
end

% Front Spar Position
fspar = [];
for i = 1:length(FAME.WING.DETAILED_GEOMETRY.FRONT_SPAR_POS.ETA)
    
    check_var = strtrim(regexp(FAME.WING.DETAILED_GEOMETRY.FRONT_SPAR_POS.ETA{i},'(?<=&).*()','match'));
    
    Eta = str2double(FAME.WING.DETAILED_GEOMETRY.FRONT_SPAR_POS.ETA{i});
    
    % If the input is defined by a user variable then call it
    if ~isempty(check_var)
        eval(['Eta = ' check_var{1} ';']);
    end
    
    fspar = [fspar; [Eta,...
        str2double(FAME.WING.DETAILED_GEOMETRY.FRONT_SPAR_POS.CHORD_REL{i})]];
    
end

% Rear Spar Position
aspar = [];
for i = 1:length(FAME.WING.DETAILED_GEOMETRY.REAR_SPAR_POS.ETA)
    
    check_var = strtrim(regexp(FAME.WING.DETAILED_GEOMETRY.REAR_SPAR_POS.ETA{i},'(?<=&).*()','match'));
    
    Eta = str2double(FAME.WING.DETAILED_GEOMETRY.REAR_SPAR_POS.ETA{i});
    
    % If the input is defined by a user variable then call it
    if ~isempty(check_var)
        eval(['Eta = ' check_var{1} ';']);
    end
    
    aspar = [aspar; [Eta,...
        str2double(FAME.WING.DETAILED_GEOMETRY.REAR_SPAR_POS.CHORD_REL{i})]];
    
end
%% Determine the number of sections for the wing
Eta = unique([unique(chord(:,1));unique(QCsweep(:,1));unique(dihedral(:,1))]);

NumSec = length(Eta)-1;

% This really should be interpolated later on!!!!
chord_eta     = interp1(chord(:,1),chord(:,2),Eta,'linear');
thickness_eta = interp1(thickness(:,1),thickness(:,2),Eta,'linear');
dihedral_eta  = interp1(dihedral(:,1),dihedral(:,2),Eta,'linear');
QCsweep_eta   = interp1(QCsweep(:,1),QCsweep(:,2),Eta,'linear');

% WING PLANFORM PROPERTIES 
param.Wing.Parent               = 'Fuselage';
param.Wing.LiftingSurface       = 1;
param.Wing.Optimisation         = 1;
param.Wing.Planform.SIMXZ       = 0;
param.Wing.Planform.NWingSec    = NumSec;
param.Wing.Planform.NBeams      = 33;

param.Wing.Planform.y           = Eta;
param.Wing.Planform.taper       = chord_eta(2:end,1)./chord_eta(1:end-1,1);

%% Interpolate the data so that it matches the kinks of the wing
param.Wing.Planform.QCsweep(:,1) = QCsweep_eta(1:end-1,1);
param.Wing.Planform.QCsweep(:,2) = QCsweep_eta(1:end-1,1);

param.Wing.Planform.Beamsweep = [];

param.Wing.Planform.LEsweep = [];

param.Wing.Planform.height_fraction(:,1) = thickness_eta(:,1);

param.Wing.Planform.dihedral(:,1) = dihedral_eta(1:end-1,1);
param.Wing.Planform.dihedral(:,2) = dihedral_eta(1:end-1,1);

param.Wing.Planform.twist(:,1) = 3*ones(size(QCsweep_eta(:,1),1),1);
param.Wing.Planform.twist(:,2) = 3*ones(size(QCsweep_eta(:,1),1),1);

param.Wing.Planform.chord       = [];
param.Wing.Planform.NumCPanels  = 1;

param.Wing.Planform.AR      = [];
param.Wing.Planform.b_ref   = str2double(FAME.WING.GLOBAL_GEOMETRY.SPAN{1})/2000;
param.Wing.Planform.c_ref   = chord(1,2);
param.Wing.Planform.S_ref   = 65; % THIS SHOULD NOT BE FIXED!!!

param.Wing.Planform.startID = 710000;

% TODO - FIX THIS SO THAT IT TAKES A DISTRIBUTION
param.Wing.Planform.fspar   = fspar(1,2);
param.Wing.Planform.aspar   = aspar(1,2);

param.Wing.Planform.t_skin0 = 0.03;
param.Wing.Planform.t_spar0 = 0.03;
param.Wing.Planform.A_str0 = 0.0001;
param.Wing.Planform.t_str0 = 0.001;

FAME_Offset = str2double(FAME.AIRPLANE.GLOBAL_DATA.WING_POSITION.X_REF{1})/1000;

Beam_Offset = 0.5*(param.Wing.Planform.aspar + param.Wing.Planform.fspar)*param.Wing.Planform.c_ref + FAME_Offset;

param.Wing.Planform.Offset = [Beam_Offset,...
    0,str2double(FAME.AIRPLANE.GLOBAL_DATA.WING_POSITION.Z_REF{1})/1000];

if isfield(FAME.WING.AERODYNAMIC_DESIGN_CASE,'LIFT_DISTRIBUTION')
    for i = 1:numel(FAME.WING.AERODYNAMIC_DESIGN_CASE.LIFT_DISTRIBUTION.ETA)
        param.Wing.Planform.Cl(i,1) = str2double(FAME.WING.AERODYNAMIC_DESIGN_CASE.LIFT_DISTRIBUTION.ETA{i});
        param.Wing.Planform.Cl(i,2) = str2double(FAME.WING.AERODYNAMIC_DESIGN_CASE.LIFT_DISTRIBUTION.CL{i});
    end
end


% CONTROL SURFACE PROPERTIES

%% Control Surface Information

param.Wing.Planform.CS.exist = 1;
count = 0;
splcount = 0;
% SPOILERS
spoiler.eta_beg_le = [];
spoiler.eta_beg_te = [];
spoiler.eta_end_le = [];
spoiler.eta_end_te = [];
spoiler.chrd_beg_le = [];
spoiler.chrd_beg_te = [];
spoiler.chrd_end_le = [];
spoiler.chrd_end_te = [];

for i = 1:length(FAME.WING.DETAILED_GEOMETRY.SPOILER_GEO.NO)
    
    count = count + 1;
    splcount = splcount + 1;
    Eta_beg_le  = str2double(FAME.WING.DETAILED_GEOMETRY.SPOILER_GEO.ETA_BEG_LE{i});
    Eta_beg_te  = str2double(FAME.WING.DETAILED_GEOMETRY.SPOILER_GEO.ETA_BEG_TE{i});
    Eta_end_le  = str2double(FAME.WING.DETAILED_GEOMETRY.SPOILER_GEO.ETA_END_LE{i});
    Eta_end_te  = str2double(FAME.WING.DETAILED_GEOMETRY.SPOILER_GEO.ETA_END_TE{i});
    Chrd_beg_le = str2double(FAME.WING.DETAILED_GEOMETRY.SPOILER_GEO.CHRD_BEG_LE{i});
    Chrd_beg_te = str2double(FAME.WING.DETAILED_GEOMETRY.SPOILER_GEO.CHRD_BEG_TE{i});
    Chrd_end_le = str2double(FAME.WING.DETAILED_GEOMETRY.SPOILER_GEO.CHRD_END_LE{i});
    Chrd_end_te = str2double(FAME.WING.DETAILED_GEOMETRY.SPOILER_GEO.CHRD_END_TE{i});
    
    spoiler.eta_beg_le = [spoiler.eta_beg_le;Eta_beg_le];
    spoiler.eta_beg_te = [spoiler.eta_beg_te;Eta_beg_te];
    spoiler.eta_end_le = [spoiler.eta_end_le;Eta_end_le];
    spoiler.eta_end_te = [spoiler.eta_end_te;Eta_end_te];
    spoiler.chrd_beg_le = [spoiler.chrd_beg_le;Chrd_beg_le];
    spoiler.chrd_beg_te = [spoiler.chrd_beg_te;Chrd_beg_te];
    spoiler.chrd_end_le = [spoiler.chrd_end_le;Chrd_end_le];
    spoiler.chrd_end_te = [spoiler.chrd_end_te;Chrd_end_te];
    
    param.Wing.Planform.CS.inboard(count) = Eta_beg_le;
    param.Wing.Planform.CS.outboard(count) = Eta_end_le;
    param.Wing.Planform.CS.chord(count) = 1 - Chrd_beg_le;
    param.Wing.Planform.CS.name{count} = ['spl' num2str(count)];
    param.Wing.Planform.CS.MaxDef(count) = 30;
    param.Wing.Planform.CS.RateLim(count) = 80;
    
end

% AILERONS
aileron.eta_beg_le = [];
aileron.eta_beg_te = [];
aileron.eta_end_le = [];
aileron.eta_end_te = [];
aileron.chrd_beg_le = [];
aileron.chrd_end_le = [];

ailcount = 0;

for i = 1:length(FAME.WING.DETAILED_GEOMETRY.AILERON_GEO.NO)
        
    count       = count + 1;
    ailcount    = ailcount + 1;
    
    Eta_beg_le  = str2double(FAME.WING.DETAILED_GEOMETRY.AILERON_GEO.ETA_BEG_LE{i});
    Eta_beg_te  = str2double(FAME.WING.DETAILED_GEOMETRY.AILERON_GEO.ETA_BEG_TE{i});
    Eta_end_le  = str2double(FAME.WING.DETAILED_GEOMETRY.AILERON_GEO.ETA_END_LE{i});
    Eta_end_te  = str2double(FAME.WING.DETAILED_GEOMETRY.AILERON_GEO.ETA_END_TE{i});
    Chrd_beg_le = str2double(FAME.WING.DETAILED_GEOMETRY.AILERON_GEO.CHRD_BEG_LE{i});
    Chrd_end_le = str2double(FAME.WING.DETAILED_GEOMETRY.AILERON_GEO.CHRD_END_LE{i});
    
    aileron.eta_beg_le = [aileron.eta_beg_le;Eta_beg_le];
    aileron.eta_beg_te = [aileron.eta_beg_te;Eta_beg_te];
    aileron.eta_end_le = [aileron.eta_end_le;Eta_end_le];
    aileron.eta_end_te = [aileron.eta_end_te;Eta_end_te];
    aileron.chrd_beg_le = [aileron.chrd_beg_le;Chrd_beg_le];
    aileron.chrd_end_le = [aileron.chrd_end_le;Chrd_end_le];
    
    param.Wing.Planform.CS.inboard(count) = Eta_beg_le;
    param.Wing.Planform.CS.outboard(count) = Eta_end_le;
    param.Wing.Planform.CS.chord(count) = 1 - Chrd_beg_le;
    param.Wing.Planform.CS.name{count} = ['ail' num2str(ailcount)];
    param.Wing.Planform.CS.MaxDef(count) = 30;
    param.Wing.Planform.CS.RateLim(count) = 80;
    
end

% MATERIAL PROPERTIES 6061-T6
% TODO - READ THIS FROM SOMEWHERE
% param.Wing.Box.E = 6.89*10^10; % Young's Modulus
% param.Wing.Box.smax = 2.76*10^8; % Bending Stress 2.76*10^8 4.75
% param.Wing.Box.shmax = 2.07*10^8; % Shear stress 2.07
% param.Wing.Box.Rho = 2.9*10^3; % Material Density
% param.Wing.Box.nu = 0.33; % Poisson's ratio
% param.Wing.Box.emax = 10^6*param.Wing.Box.smax/param.Wing.Box.E;
% param.Wing.Box.smax = 2.76*10^8; % Bending Stress 2.76*10^8
% param.Wing.Box.shmax = 0.55*2.76*10^8; % Bending Stress 2.76*10^8
param.Wing.Box.KAs = 0.36;
param.Wing.Box.ds = 0.12; % Stringer thickness/depth
param.Wing.Box.SP = 0.15; % Stringer Pitch

param.Wing.Box.E = 7.8437*10^10; % Young's Modulus 7.843720000000000e
param.Wing.Box.emax = 2*3500/1.5; 
param.Wing.Box.smax = param.Wing.Box.emax*param.Wing.Box.E/(10^6); % Bending Stress 2.76*10^8
param.Wing.Box.shmax = 1.8019e+08; % Shear stress 4*10^7
param.Wing.Box.Rho = 1.6*10^3; % Material Density
param.Wing.Box.nu = 0.4978; % Poisson's ratio

if ~isfield(FAME.USER_VARIABLES,'rib_spc')
    fprintf('\nWARNING - NO RIB SPACING FOUND. ASSIGNED VALUE OF RP = 0.65m');
    param.Wing.Box.RP = 0.65;
else
    param.Wing.Box.RP = str2double(FAME.USER_VARIABLES.rib_spc{1,1})/1000; %
end

% CARRYTHROUGH SECTION
param.Wing.Planform.CT_y_eta = root_wing;

%% VTP PROPERTIES

param.VTP.LiftingSurface = 1;
param.VTP.Parent = 'Fuselage';
param.VTP.Optimisation   = 0;
param.VTP.Planform.SIMXZ = 1;

param.VTP.Planform.NWingSec = 1;
param.VTP.Planform.NBeams = 5;

param.VTP.Planform.y = [0.0;1.0];
param.VTP.Planform.taper = [0.9985];

param.VTP.Planform.QCsweep(:,1) = 37.4879;
param.VTP.Planform.QCsweep(:,2) = 37.4879;

param.VTP.Planform.Beamsweep = [];

param.VTP.Planform.LEsweep = [];

param.VTP.Planform.height_fraction = [0.03;0.03];

param.VTP.Planform.dihedral(:,1) = 90;
param.VTP.Planform.dihedral(:,2) = 90;

param.VTP.Planform.twist(:,1) = 0;
param.VTP.Planform.twist(:,2) = 0;

param.VTP.Planform.chord = [];
param.VTP.Planform.NumCPanels = 1;

param.VTP.Planform.AR    = [];
param.VTP.Planform.b_ref = 5.875;
param.VTP.Planform.c_ref = 5.0354;
param.VTP.Planform.S_ref = [];

param.VTP.Planform.Offset = [38.76,0.0,2.1145] + [0.5*(0.1+0.65)* 5.0354,0,0];
param.VTP.Planform.startID = 910000;
param.VTP.Planform.fspar = 0.10;
param.VTP.Planform.aspar = 0.65;

% CONTROL SURFACE PROPERTIES
param.VTP.Planform.CS.exist = 1;
param.VTP.Planform.CS.inboard = 0;
param.VTP.Planform.CS.outboard = 0.9;
param.VTP.Planform.CS.chord = 0.35;
param.VTP.Planform.CS.name{1} = 'rudr';

% MATERIAL PROPERTIES
param.VTP.Box.E = 6.89*10^15; % Young's Modulus
param.VTP.Box.smax = 2.80*10^8; % Bending Stress
param.VTP.Box.shmax = 2.80*10^8; % Shear stress
param.VTP.Box.Rho = 2.9*10^3; % Material Density
param.VTP.Box.nu = 0.33; % Poisson's ratio
param.VTP.Box.emax = 10^6*param.Wing.Box.smax/param.Wing.Box.E; 

param.VTP.Box.KAs = 0.36;
param.VTP.Box.ds = 0.12; % Stringer depth
param.VTP.Box.SP = 0.15; % Stringer Pitch
param.VTP.Box.RP = 0.6; % Rib Pitch

% ENGINE PROPERTIES
param.VTP.Engine.Exist = 1;
param.VTP.Engine.y = 0.0;
param.VTP.Engine.offset = [40.57,3.9,4.1117] - [38.3186,0.0,1.8198] - [2.7880,0,0];
param.VTP.Engine.mass = 3440.58;
param.VTP.Engine.Planform.Length = 5.57604;
param.VTP.Engine.Planform.y_eta = [0;0.5;0.8;1];
param.VTP.Engine.Planform.radius = 3.044/2*[1;1;0.9;0.8];
param.VTP.Engine.Planform.NSec = length(param.VTP.Engine.Planform.y_eta)-1;

param.VTP.Planform.t_skin0 = 0.01;
param.VTP.Planform.t_spar0 = 0.01;
param.VTP.Planform.A_str0 = 0.0001;
param.VTP.Planform.t_str0 = 0.001;

%% HTP PLANFORM PROPERTIES

if ~isfield(FAME.AIRPLANE,'HORIZONTAL_TAILPLANE_POSITION')
    
    param.HTP.Planform.NBeams = 5;
    
    %% EVALUATE USER VARIABLES
    
    %% Extract some key wing properties from the FAME file
    
    
    % Quarter Chord Sweep
    QCsweep = [];
    for i = 1:length(FAME.HORIZONTAL_TAILPLANE.GLOBAL_GEOMETRY.SWEEP_DISTRIBUTION.SWEEP_25)
        
        check_var = strtrim(regexp(FAME.HORIZONTAL_TAILPLANE.GLOBAL_GEOMETRY.SWEEP_DISTRIBUTION.ETA_BEG{i},'(?<=&).*()','match'));
        
        Eta = str2double(FAME.HORIZONTAL_TAILPLANE.GLOBAL_GEOMETRY.SWEEP_DISTRIBUTION.ETA_BEG{i});
        % If the input is defined by a user variable then call it
        if ~isempty(check_var)
            eval(['Eta = ' check_var{1} ';']);
        end
        
        QCsweep = [QCsweep; [Eta,...
            str2double(FAME.HORIZONTAL_TAILPLANE.GLOBAL_GEOMETRY.SWEEP_DISTRIBUTION.SWEEP_25{i})]];
        
    end
    QCsweep(end+1,:) = [1,QCsweep(end,2)];
    
    % Thickness to Chord Ratio
    thickness = [];
    for i = 1:length(FAME.HORIZONTAL_TAILPLANE.GLOBAL_GEOMETRY.THICKNESS_DISTRIBUTION.THICKNESS_REL)
        
        check_var = strtrim(regexp(FAME.HORIZONTAL_TAILPLANE.GLOBAL_GEOMETRY.THICKNESS_DISTRIBUTION.ETA{i},'(?<=&).*()','match'));
        
        Eta = str2double(FAME.HORIZONTAL_TAILPLANE.GLOBAL_GEOMETRY.THICKNESS_DISTRIBUTION.ETA{i});
        
        % If the input is defined by a user variable then call it
        if ~isempty(check_var)
            eval(['Eta = ' check_var{1} ';']);
        end
        
        thickness = [thickness; [Eta,...
            str2double(FAME.HORIZONTAL_TAILPLANE.GLOBAL_GEOMETRY.THICKNESS_DISTRIBUTION.THICKNESS_REL{i})]];
        
    end
    
    % Dihedral
    dihedral = [];
    for i = 1:length(FAME.HORIZONTAL_TAILPLANE.GLOBAL_GEOMETRY.DIHEDRAL_DISTRIBUTION.DIHEDRAL)
        
        check_var = strtrim(regexp(FAME.HORIZONTAL_TAILPLANE.GLOBAL_GEOMETRY.DIHEDRAL_DISTRIBUTION.ETA_BEG{i},'(?<=&).*()','match'));
        
        Eta = str2double(FAME.HORIZONTAL_TAILPLANE.GLOBAL_GEOMETRY.CHORD_DISTRIBUTION.ETA{i});
        
        % If the input is defined by a user variable then call it
        if ~isempty(check_var)
            eval(['Eta = ' check_var{1} ';']);
        end
        
        dihedral = [dihedral; [Eta,...
            str2double(FAME.HORIZONTAL_TAILPLANE.GLOBAL_GEOMETRY.DIHEDRAL_DISTRIBUTION.DIHEDRAL{i})]];
        
    end
    dihedral(end+1,:) = [1,dihedral(end,2)];
    
    % Chord Distribution
    chord = [];
    for i = 1:length(FAME.HORIZONTAL_TAILPLANE.GLOBAL_GEOMETRY.CHORD_DISTRIBUTION.CHORD)
        
        check_var = strtrim(regexp(FAME.HORIZONTAL_TAILPLANE.GLOBAL_GEOMETRY.CHORD_DISTRIBUTION.ETA{i},'(?<=&).*()','match'));
        
        Eta = str2double(FAME.HORIZONTAL_TAILPLANE.GLOBAL_GEOMETRY.CHORD_DISTRIBUTION.ETA{i});
        
        % If the input is defined by a user variable then call it
        if ~isempty(check_var)
            eval(['Eta = ' check_var{1} ';']);
        end
        
        chord = [chord; [Eta,...
            str2double(FAME.HORIZONTAL_TAILPLANE.GLOBAL_GEOMETRY.CHORD_DISTRIBUTION.CHORD{i})/1000]];
        
    end
    
    %% Determine the number of sections for the wing
    Eta = unique([unique(chord(:,1));unique(QCsweep(:,1));unique(dihedral(:,1))]);
    
    NumSec = length(Eta)-1;
    
    param.HTP.Planform.y        = Eta;
    param.HTP.Planform.NWingSec = NumSec;
    
    % This really should be interpolated later on!!!!
    chord_eta     = interp1(chord(:,1),chord(:,2),Eta,'linear');
    thickness_eta = interp1(thickness(:,1),thickness(:,2),Eta,'linear');
    dihedral_eta  = interp1(dihedral(:,1),dihedral(:,2),Eta,'linear');
    QCsweep_eta   = interp1(QCsweep(:,1),QCsweep(:,2),Eta,'linear');
    
    param.HTP.Planform.taper = chord_eta(2:end,1)./chord_eta(1:end-1,1);
    param.HTP.LiftingSurface  = 1;
    param.HTP.Parent  = 'Fuselage';
    param.HTP.Planform.SIMXZ       = 0;
    param.HTP.Optimisation = 0;
    
    %% Interpolate the data so that it matches the kinks of the wing
    param.HTP.Planform.height_fraction(:,1) = thickness_eta(:,1);
    
    param.HTP.Planform.dihedral(:,1) = dihedral_eta(1:end-1,1);
    param.HTP.Planform.dihedral(:,2) = dihedral_eta(1:end-1,1);
    
    param.HTP.Planform.QCsweep(:,1) = QCsweep_eta(1:end-1,1);
    param.HTP.Planform.QCsweep(:,2) = QCsweep_eta(1:end-1,1);
    
    param.HTP.Planform.Beamsweep = [];
    
    param.HTP.Planform.LEsweep = [];
    
    param.HTP.Planform.twist(:,1) = 3*ones(size(QCsweep_eta(:,1),1),1);
    param.HTP.Planform.twist(:,2) = 3*ones(size(QCsweep_eta(:,1),1),1);
    
    param.HTP.Planform.chord = [];
    param.HTP.Planform.NumCPanels = 1;
    %
    param.HTP.Planform.AR = [];
    param.HTP.Planform.b_ref = str2double(FAME.HORIZONTAL_TAILPLANE.GLOBAL_GEOMETRY.SPAN{1})/2000;
    param.HTP.Planform.c_ref = chord(1,2);
    
    param.HTP.Planform.S_ref = [];
    
    param.HTP.Planform.Offset = [42.884,0,7.061] + [0.5*(0.15+0.65)*chord(1,2),0,0]; % [16,0,0]
    
    param.HTP.Planform.startID = 810000;
    param.HTP.Planform.fspar = 0.15;
    param.HTP.Planform.aspar = 0.65;
    
    param.HTP.Planform.t_skin0 = 0.01;
    param.HTP.Planform.t_spar0 = 0.01;
    param.HTP.Planform.A_str0 = 0.0001;
    param.HTP.Planform.t_str0 = 0.001;
    
    % CONTROL SURFACE PROPERTIES
    param.HTP.Planform.CS.exist = 1;
    param.HTP.Planform.CS.inboard(1) = 0.05;
    param.HTP.Planform.CS.outboard(1) = 1;
    param.HTP.Planform.CS.chord(1) = 0.35;
    param.HTP.Planform.CS.name{1} = 'elv';
    param.HTP.Planform.CS.MaxDef(1) = 30;
    param.HTP.Planform.CS.RateLim(1) = 80;
    
    % MATERIAL PROPERTIES 6061-T6
    param.HTP.Box.E = 6.89*10^15; % Young's Modulus
    param.HTP.Box.smax = 2.76*10^8; % Bending Stress 2.76*10^8 4.75
    param.HTP.Box.shmax = 2.07*10^8; % Shear stress 2.07
    param.HTP.Box.Rho = 1.1*10^3; % Material Density
    param.HTP.Box.nu = 0.33; % Poisson's ratio
    param.HTP.Box.emax = 10^6*param.Wing.Box.smax/param.Wing.Box.E;
    % param.Wing.Box.smax = 2.76*10^8; % Bending Stress 2.76*10^8
    % param.Wing.Box.shmax = 0.55*2.76*10^8; % Bending Stress 2.76*10^8
    param.HTP.Box.KAs = 0.36;
    param.HTP.Box.ds = 0.12; % Stringer thickness/depth
    param.HTP.Box.SP = 0.15; % Stringer Pitch
    
    if ~isfield(FAME.USER_VARIABLES,'rib_spc')
        fprintf('\nWARNING - NO RIB SPACING FOUND. ASSIGNED VALUE OF RP = 0.65m');
        param.HTP.Box.RP = 0.65; %
    else
        param.HTP.Box.RP = str2double(FAME.USER_VARIABLES.rib_spc{1,1})/1000; %
    end

else
    param.HTP = [];
end

%% FUSELAGE PROPERTIES
param.Fuselage.LiftingSurface  = 0;
param.Fuselage.Parent  = [];

param.Fuselage.Planform.NSec = 40;

param.Fuselage.Planform.Length = 43.913;
fprintf('\n Warning ... The Fuselage length has been defined as being 44 m .');
fprintf('\n This will have to be modified so that the fuselage length is extracted');
fprintf('\n from the some FAME input file.');

param.Fuselage.Planform.y_eta = [0;0.005268;0.010536;0.015805;0.021073;...
    0.026342;0.03161;0.036879;0.042147;0.047415;0.052684;...
    0.057952;0.063221;0.0684890;0.073758;0.079026;0.084294;0.089563;0.094831;0.1001;0.411022;...
    0.721944;0.736578;0.751213;0.765847;0.780482;0.795117;0.809751;0.824386;0.83902;0.853655;...
    0.868289;0.882924;0.897558;0.912193;0.926827;0.941462;0.956096;0.970731;0.985365;1];

param.Fuselage.Planform.radius = [0;0.3844030;0.565081;0.707928;0.830682;0.940375;...
    1.04067;1.13377;1.22112;1.30374;1.38237;1.45758;1.52981;1.59941;1.66667;...
    1.73182;1.79508;1.8566;1.91653;1.975;2.11455;2.11455;2.11455;2.11455;2.11455;...
    2.11455;2.11455;2.11455;2.11455;2.11455;2.11455;2.11455;2.11455;2.11455;1.9;...
   1.75;1.6;1.4;1.2;1.0;0];

param.Fuselage.Planform.Offset = [0,0,0];

param.Fuselage.Planform.Mass = 24126;
param.Fuselage.Planform.CG   = [22,0,0];

% MATERIAL PROPERTIES
param.Fuselage.Box.E = 6.89*10^10; % Young's Modulus
param.Fuselage.Box.smax = 2.80*10^8; % Bending Stress
param.Fuselage.Box.shmax = 2.80*10^8; % Shear stress
param.Fuselage.Box.Rho = 2.9*10^3; % Material Density
param.Fuselage.Box.nu = 0.33; % Poisson's ratio
param.Fuselage.Box.emax = 10^6*param.Wing.Box.smax/param.Wing.Box.E;

param.Fuselage.Box.KAs = 0.36;
param.Fuselage.Box.ds = 0.12; % Stringer depth
param.Fuselage.Box.SP = 0.15; % Stringer Pitch
param.Fuselage.Box.RP = 0.6; % Rib Pitch
param.Fuselage.Properties.Mass = 24126;
param.Fuselage.Properties.CG   = [22,0,0];

% PAYLOAD PROPERTIES
param.Payload.LiftingSurface = 0;
param.Payload.Parent = 'Fuselage';

param.Payload.Planform.y_eta  = [0;1];
param.Payload.Planform.radius = [1.5;1.5];
param.Payload.Planform.Offset = [6.4,0,0];
%param.Payload.Planform.GlobalOffset = [16.4,0,0];
param.Payload.Planform.Length = 20;
param.Payload.Planform.NSec    = 1;
%param.Payload.Properties.Mass   = 16059/2;
param.Payload.Properties.Mass   = 2.92e+04 / 2;
param.Payload.Properties.CG     = [16.4,0,0];
param.Payload.Planform.Colour = [0,1,0];
    
% ENGINE PROPERTIES
param.Engine.LiftingSurface = 0;
param.Engine.Parent = 'VTP';

param.Engine.Planform.NSec = 3;
param.Engine.Planform.Length = 5.57604;
param.Engine.Planform.y_eta  = [0;0.5;0.8;1];
param.Engine.Planform.radius = 3.044*[1;1;0.9;0.8]/2;
param.Engine.Planform.Offset = [36.7,3.9,4.1117] - [40.6483,0.0,0.0] - [5.57604/2,0,0];
%param.Engine.Planform.GlobalOffset = [36.7,3.9,4.1117];
param.Engine.Planform.y_eta_parent = 0.0;

param.Engine.Properties.Mass = (9525+776)/2;
param.Engine.Properties.TSFC = 0.5649/(60*60);
param.Engine.Properties.CG     = [36.7,3.9,4.1117];
param.Engine.Planform.Colour = [1,0,0];

check_var = strtrim(regexp(FAME.AIRPLANE.GLOBAL_DATA.DESIGN_WEIGHTS.MTOW{1},'(?<=&).*()','match'));

% If the input is defined by a user variable then call it
if ~isempty(check_var)
    eval(['MTOW = ' check_var{1} ';']);
else
    MTOW = str2double(FAME.AIRPLANE.GLOBAL_DATA.DESIGN_WEIGHTS.MTOW{1});
end

%% WEIGHT PROPERTIES
param.Weights.MTOW       = MTOW;
param.Weights.OEW        = 51967;
param.Weights.FB         = 14800;
param.Weights.Payload    = 2.92e+04;
end