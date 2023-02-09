%% Wing model
Wing = awi.model.LiftingSurface;

% define number of beam element
Wing.NumBeamElem = 10;

% set model name
Wing.Name = 'Wing_model';

% set model origin - root
Wing.Origin=[0,0,0]; 

%Use the Leading/Trailing edge sweep to define the planform
Wing.ActiveSet = 'sSet';

% Wing planform dimensions
Wing.SpanVector  = 'Y';
Wing.Span        = 10;

Wing.LESweep     = [0,0];
Wing.LESweep_eta = [0,1];

Wing.TESweep     = [0,0];
Wing.TESweep_eta = [0,1];

Wing.RootChord   = 1;

% Dihedral
Wing.Dihedral=[0,0];
Wing.Dihedral_eta=[0,1];

% Make sure the beam is at the midchord hence 0.5
all_eta           = Wing.Eta_;
Wing.BeamLoc     = repmat(0.5, size(all_eta));
Wing.BeamLoc_eta = all_eta;

%Make the material
E_wing  = 71.7e9; %[N/m^2], typical YM of aluminium
nu_wing = 0.333;
rho_wing=2810;

Mat_wing = awi.model.Material;
Mat_wing.E  = E_wing;
Mat_wing.Nu = nu_wing;
Mat_wing.G  = E_wing / (2 * (1 + nu_wing));
Mat_wing.Rho=rho_wing;

Wing.Material_eta = [0, 1];
Wing.Material     = [Mat_wing, Mat_wing];

% AeroPanelLength
Wing.AeroPanelLength=[];
Wing.NumAeroPanel=10;
Wing.AeroPanelAR=1;

%% set a genetric material properties 
Wing.A   =  [0.005,0.005];     % values 
Wing.A_eta = [0,1];            % relative nodes position 0 at root; 1 at tip.

Wing.I11 = [0.1,0.1];
Wing.I11_eta = [0,1];

Wing.I22 = [0.01,0.01];
Wing.I22_eta = [0,1];

Wing.J   = [0.01,0.01];
Wing.J_eta= [0,1];

build(Wing)

% convert to FEM 
FEM_wing = convertToFE(Wing);
draw(FEM_wing)

% write to file
if ~isfolder(fullfile(pwd,'bin'))
   mkdir(fullfile(pwd,'bin'))
else
   delete(fullfile('bin','*'))
end

folder = fullfile(pwd,'bin'); 

[~, includeFiles] = export(FEM_wing, folder, 'WriteHeaderFile', false);
