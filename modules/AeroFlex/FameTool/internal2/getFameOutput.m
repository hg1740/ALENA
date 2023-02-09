%%GETFAMEOUTPUT:
% The following function can take an object from the FAME converstion tool
% and populae it with both the results and the filenames associated to
% them. If in the case that no object is passed into the function,
% getFameOutput queries the user for the FAME folder path.
%
%   Author: Dario Calderon (University of Bristol 2017)

function obj = getFameOutput(obj, logfcn)

if nargin < 2
    logfcn = @(s) fprintf('%s\n', s);
end

% Conditional statement in the case that there is no input object
if nargin < 1
    DirName = uigetdir('Locate the master FAME folder for the model');
    obj.Dirs.fameFolder = DirName;
    obj.Inp = [];
    objType = 0;
else
    objType = 1;
end

% Check if the export folder exists
exportexists = isdir([obj.Dirs.fameFolder,filesep,'export']);

if exportexists == 1
    
    % Read INTERNAL LOADS results
    [Results,inputLoadFiles]        = getFameInternalLoads(obj.Dirs.fameFolder);
    obj.Fame.Results.InternalLoads  = Results;
    
    % Read AERODYNAMIC results
    [Results,inputAeroFiles]        = getFameAerodynamics(obj.Dirs.fameFolder);
    obj.Fame.Results.Aerodynamic    = Results;
    
    % Read BOX PROPERTIES results
    [Results,inputThickFiles]       = getFameThickness(obj.Dirs.fameFolder);
    obj.Fame.Results.Thickness      = Results;
    
    % Read DEFORMATION results - The F
    [Results,inputDefFiles]         = getFameDeformations(obj.Dirs.fameFolder);
    obj.Fame.Results.Deformation    = Results;
    
    % Read BOX STIFFNESS properties
    [Results,inputStiffFiles]       = getFameStiffness(obj.Dirs.fameFolder);
    obj.Fame.Results.Stiffness      = Results;
    
    % Read BOX AREAS
    [Results,inputBoxFiles]         = getFameBoxAreasFcn(obj.Dirs.fameFolder);
    obj.Fame.Results.BoxAreas       = Results;
    
    % % Read Stresses
    % [Results,inputStressFiles]      = getFameStress(obj.Dirs.fameFolder);
    % obj.Fame.Results.Stress         = Results;
     
else
    inputLoadFiles.internalLoadsFiles{1}    = '';
    inputAeroFiles.aeroCasesFiles{1}        = '';
    inputThickFiles.thicknessFiles{1}       = '';
    inputDefFiles.wingDeformationsFiles{1}  = '';
    inputStiffFiles.boxPropertiesFiles{1}   = '';
    inputBoxFiles.boxAreaFiles{1}           = '';
    
    logfcn(' * * WARNING * * - Fame Outputs not found');
    
end

% Read LOADCASE

% Check to see if the 1_structure file exists, if it doesn't 
lcexists = isdir([obj.Dirs.fameFolder,filesep,'1_structure']);

if lcexists == 1
    
    pathName = ['1_structure',filesep,'1_10_wing',filesep,'l3_fame-w',filesep,'flexible'];
    fileName = 'loadcase_info.txt';
    
    filelist = dir(fullfile(obj.Dirs.fameFolder,pathName,fileName));
    
    if isempty(filelist)
        [fName,pName,idx] = uigetfile([pwd '*txt'],'Pick the loadcase_info.txt file');
        if idx~= 0
            obj.Inp.loadCaseFiles{1} = fullfile(pName,fName);
        else
            error('Unable to proceed without the loadcase_info.txt file');
        end
    else
        obj.Inp.loadCaseFiles{1} = fullfile(obj.Dirs.fameFolder,pathName,fileName);
    end
    
    [Results,inputLCFiles] = getFameLoadCases(obj.Inp.loadCaseFiles{1});
    
    obj.Fame.LoadCases     = Results;
else
    obj.Inp.loadCaseFiles{1} = '';
    obj.Fame.LoadCases     = [];
end


% Pass the file paths and names back to the object/structure
if objType == 1
    obj.Inp = setFields(obj.Inp,inputLoadFiles);
    obj.Inp = setFields(obj.Inp,inputAeroFiles);
    obj.Inp = setFields(obj.Inp,inputThickFiles);
    obj.Inp = setFields(obj.Inp,inputDefFiles);
    obj.Inp = setFields(obj.Inp,inputStiffFiles);
    obj.Inp = setFields(obj.Inp,inputBoxFiles);
    %obj.Inp = setFields(obj.Inp,inputStressFiles);
else
    %obj = setfield(obj,'Inp');
end

end