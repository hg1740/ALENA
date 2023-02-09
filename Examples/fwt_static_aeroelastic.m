%% fwt_static_aeroelastic
% 
% Runs a NASTRAN static aeroelastic analysis of a cantilevered wing model
% with a folding wing tip.
%
%
% Next steps:
%   - Set the `HingeStiffness` of the FWT to a sensible number
%   - Run a normal modes analysis of the FWT
%   - Run a static aeroelastic analysis (fixed AoA) for the wing only
%   - Run a static aeroelastic analysis (fixed AoA) for FWT
%   - Set the `HingeStiffness` of the FWT to a sensible number

fold_angle  = 0;   %[deg],
flare_angle = -30;   %[deg],
aoa = 5;

if ~isfolder(fullfile(pwd,'bin'))
   mkdir(fullfile(pwd,'bin'))
else
   delete(fullfile('bin','*'))
end
run_folder = fullfile(pwd,'bin');    

acMass = 500;

altitude          = 0;
mach_number       = 0.3;
aircraft_velocity = 40;

%% Generate the geometry objects

% Make the model and visualise it
%   - Use default Hodges wing
LS = awi.model.LiftingSurface.makeHodgesWing;
LS.Twist = deg2rad([aoa,aoa]);
LS.Twist_eta = [0,1];
%Build the object in order to populate remaining properties
build(LS);


% Add a wing fold at the desired flare angle 
%   - Angles in degrees
FWT = insertWingFold(LS,'EtaFold',0.75, 'FlareAngle', flare_angle, 'FoldAngle', fold_angle);
FWT.HingeStiffness = [0 0 0 1e-4 0 0];
FWT.Twist = repmat(atan(sind(flare_angle)*sind(fold_angle)),1,2)+cosd(fold_angle)*deg2rad(aoa);
FWT.Twist = deg2rad([5,-15]);
FWT.Twist_eta = [0,1];
build(LS);
%The analysis methods require an 'awi.model.Aircraft' object
% This is because some information is only known at the aircraft level,
% e.g. all-up mass, reference span, area, etc.
Aircraft = awi.model.Aircraft;
Aircraft.add(LS);

%Set up reference quantities
%   - Note that as we have a wing fold we need to combine the properties of
%     the parent and child objects
Aircraft.RefArea  = sum([LS.SurfaceArea, FWT.SurfaceArea]);
Aircraft.RefSpan  = LS.Span + FWT.Span;
Aircraft.RefChord = LS.RootChord;

%% Generate the FEM

% Convert to a finite element model and draw it
FEM = convertToFE(Aircraft);
FEM.updateDmiEntry();
draw(FEM);

% sort out DMI card
val = flatlist(FEM);

% %Export it to a file
export(FEM, run_folder);

%% Run the analysis
TrimLoadcase = awi.model.LoadCase;
TrimLoadcase.Name = 'Fixed-semiwing';
TrimLoadcase.Altitude   = altitude;
TrimLoadcase.Mach       = mach_number;
TrimLoadcase.AcVelocity = aircraft_velocity;
TrimLoadcase.AcMass = 1;
TrimLoadcase.LoadFactor = 1;
% 
build(TrimLoadcase)

% Write input and run analysis

NastranMethods1 = awi.methods.Nastran;
NastranMethods1.AnalysisModel = FEM;
MassCases = [];

trimFile1 = NastranMethods1.writeTrimFile(Aircraft, TrimLoadcase,...
    MassCases,run_folder,'DatFilename','fwt','AoA',deg2rad(aoa));

NastranMethods1.runNastran(trimFile1);

model = mni.import_matran(fullfile(run_folder,'fwt.dat'));
model.draw

%extract the data
f06 = mni.result.f06(fullfile(run_folder,'fwt.f06'));
res_disp =  f06.read_disp;
res_aeroP = f06.read_aeroP;
res_aeroF = f06.read_aeroF;

% apply deformation result
[~,i] = ismember(model.GRID.GID,res_disp.GP);
model.GRID.Deformation = [res_disp.dX(:,i);res_disp.dY(:,i);res_disp.dZ(:,i)];

% apply aero pressure
model.CAERO1.PanelPressure = res_aeroP.Cp;

%apply aero forces
f = [res_aeroF.aeroFx;res_aeroF.aeroFy;res_aeroF.aeroFz;...
    res_aeroF.aeroMx;res_aeroF.aeroMy;res_aeroF.aeroMz];
model.CAERO1.PanelForce = f';

% update the plot to apply deformations and aero pressures + forces
model.update('Scale',1)

