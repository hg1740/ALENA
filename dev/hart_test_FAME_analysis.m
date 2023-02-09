function hart_test_FAME_analysis(AlenaFrmwrk)
%hart_test_FAME_analysis Runs some generic analysis of the HARTEN FAME
%model to test the functionality of ALENA.

close all

if nargin < 1 %Import HARTEN model as standard
    fm4File = ['\\ads\filestore\Engineering\Research\Projects\AWI\', ...
        'Analysis_Tools\UoB_Framework\Models\HARTEN_v2_VLM\06_weights\', ...
        'FAME\v3_newCG_LC.fame-w.4.00.m1-1-pv8\__mirror\newCG_LC.fm4'];
    AlenaFrmwrk = awi.model.Framework;
    AlenaFrmwrk = import(AlenaFrmwrk, fm4File, '-batch');
end

%% Scripting the analysis
% FEM = convertToFE(AlenaFrmwrk.Aircraft);
% 
% %Use the Intrinsic, Strain-based FE formulation
% MethodObj = awi.methods.IntrinsicStrainFE;
% MethodObj.AnalysisModel = FEM;


%% Old method

%Run Robbie's code without going through GLT package
% BeamModel = findall(AlenaFrmwrk, 'Type', 'BeamModel'); 
% [ISBNFE.Matrices, ISBNFE.Aero,~] = InitialiseRobbie(BeamModel.BM);
% [ISBNFE.Matrices, ISBNFE.Aero,~] = InitialiseRobbie_old(BeamModel.BM);

%Run Dario's code without going through GLT
% runSizeAnalysisNoGLT(AlenaFrmwrk);

%Makr the sizing view and run a sizing analysis using Dario's code
% SizeView = awi.view.SizeAnalysis(AlenaFrmwrk);
% analyse(SizeView);

%Make the trim view and run a trim analysis using Robbie's code
TrimView = awi.view.TrimAnalysis(AlenaFrmwrk);
x_trim = analyse(TrimView, 'AnalysisMethod', 'Instrinsic Strain-Based FE');

%Check the results
baseline_x_trim = load('BaselineTrimResults.mat', 'x_trim');
baseline_x_trim = baseline_x_trim.x_trim;

assert(all(arrayfun(@(i) isequal(x_trim.x_disp{i}, x_trim.x_disp{i}), ...
    1 : numel(baseline_x_trim.x_disp))), ...
    '* * * ERRROR IN THE TRIM RESULTS!!! CHECK CODE. * * *');

fprintf('\n\n\tANALYSIS SUCCESSFUL!!\n\n');

end


function runSizeAnalysisNoGLT(Alena)

nElem = 25;

%Find all lifting surfaces and bluff bodies
LSObj = findall(Alena.Aircraft, 'Type', 'LiftingSurface');
BBObj = findall(Alena.Aircraft, 'Type', 'BluffBody');

%BUT convertFmwk2Neo ASSUMES that all lifting surfaces have spars and control surfaces,
% so to avoid downstream error, identify and eliminate any that do not
%TODO: Check with Robbie what is right thing to do here
LSObj(arrayfun(@(x)isempty(findall(x, 'Type', 'Spar')), LSObj)) = [];
LSObj(arrayfun(@(x)isempty(x.ControlSurfaces), LSObj)) = [];

% Convert to a usable Aeroflex workspace
[beam_model,Aircraft_param,Optim_Variables,LoadCases] = ...
    convertFmwk2Neo(LSObj,BBObj,Alena.Aircraft,Alena.LoadCases,nElem);

% Restructure the variables so that they can be used by the
% AircraftSizer
GUI_Input                = awi.view.SizeAnalysis.InitialiseGUI;
GUI_Input.LoadCases      = LoadCases;
GUI_Input.beam_model     = beam_model;
GUI_Input.Aircraft_param = Aircraft_param;

GUI_Input.Optim_Variables = Optim_Variables;

% % Initialise the panel for the view
% hp = uiextras.TabPanel('Parent', obj.hPanel, 'TabWidth', 100);
% hb = uiextras.HBox('Parent', hp);
% hp.TabTitles{end} = 'Progress';

% Update the Aeroflex variable according to the
% linear/nonlinear choice

GUI_Input.Parameters.Struc_optimisation = true;
GUI_Input.Parameters.Aero_optimisation  = false;

GUI_Input.Parameters.AeroDistribution = 'Elliptical';

bLinear = true;

if bLinear
    
    GUI_Input.Parameters.Linear = true;
    
%     obj.AnalysisType{end+1} = 'Lin';
    
    % Run the sizing
    [TrimResults,OptimVariables,AeroFlexModel] = AircraftSizer(GUI_Input);
    
    % Pass the beam model here to extract what is necessary for
    % the framework
%     obj.BeamModel{end + 1} = ExtractBareBeamModel(obj,AeroFlexModel);

else
    
    GUI_Input.Parameters.Linear = false;
    
    obj.AnalysisType{end+1} = 'NonLin';
    
    % Run the sizing
    [obj.TrimResults{end+1},obj.OptimVariables{end+1},AeroFlexModel] = AircraftSizer(GUI_Input,hb, logfcn);
    
    % Pass the beam model here to extract what is necessary for
    % the framework
    obj.BeamModel{end + 1} = ExtractBareBeamModel(obj,AeroFlexModel);
    
end


end


