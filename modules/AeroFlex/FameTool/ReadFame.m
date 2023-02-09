%READFAME: A FAME model conversion script
%
% Summary: The following script looks through the FAME folder, reads the
% data and passes it into the matlab workspace. This can then be used used
% to generate NASTRAN/ADAMS Input file
%
% Date of Creation: 12/2015
%
% Authors:  Dario Calderon - University of Bristol
%           Etienne Coetzee - Airbus (FPO)
%
%
%% Add subfolders to the search path
set_path % TODO - I added this script to skip a step - Dario -> Etienne

%% Create conversion object first
% We are using objects to handle the data, so let us beging. Create the
% initial object
a_fame = Fame2mat;

%% Set directories
% Make sure that you directory settings are correct. You need to set the
% following
%
% * home folder  - location of the NEOCASS program
% * FAME folder  - location of FAME output files
% * Fuel data    - location of output after running FAME TANK
% * write folder - location where you want output files written
a_fame.Dirs

%% Properties
% Now you can check to see what is in this object. Type 
%
% |>> a|
%
% at the command line or |properties(a)| .
%
% You can see that there are 3 property fields in this object
% 
% * Dirs - object that contains paths to model directories etc
% * Opts - object that contains options for meshing, aerodynamics etc
% * Mdl  - model object that contains all the structural and aerodynamic
%          information, as well as load case information.
% 
% Note that there are some hidden properties that can be used by an expert
% user, but we decided to hide these. You can see what these are by looking
% in the class definition file for this object.
properties(a_fame)

%% Creat mass corrections file
% In some cases we need adjust the mass data, hence we create a mass
% corrections file. We need to tell the program that this file is present.
% We do this by setting the input flag and defining the path to the
% corrections file
a_fame.Opts.Struct.massCorrection = true;

folderName              = a_fame.Dirs.fameFolder;
wbdPathName             = ['__mon',filesep,'wbd2',filesep,'massdist',filesep,'flexible'];

a_fame.Inp.wbdCorrFiles = fullfile(folderName,wbdPathName,'mass_corrections.xlsx');

%% Methods
% Each method can be called on its own, but the FAME data should be
% extracted first. Let us see what methods are available.
methods(a_fame)

%% Airfoil Profiles
% This generates all the airfoil profiles for the the root and tip of each
% panel
%interpprofneocassformat(a_fame.Opts.Aero.nSpan);
%% Run the model
% The |run| method is used to run through the whole process.
a_fame = run(a_fame);

%% Now plot the structural mesh
plotNeocassModel(a_fame);

Fame = a_fame.Fame;

save([a_fame.Dirs.writeFolder filesep 'HartenV2'],'Fame');

clear ans Fame folderName wbdPathName;


