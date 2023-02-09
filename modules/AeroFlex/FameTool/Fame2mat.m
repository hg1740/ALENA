%FAME2MAT: Object generator for the ReadFame Script
%
% Summary:
%
% Date of Creation: 12/2015
%
% Authors:  Dario Calderon  - University of Bristol
%           Etienne Coetzee - Airbus (FPO)
%
%
classdef Fame2mat < hgsetget
    %% Properties
    % Can be accesed by the user
    properties(Access = public)
        Dirs            % Directory Options
        Opts            % General run-time Options
        Mdl             % Information for FE model based on FAME model
    end
    
    %Logical flags indicating which files/directories the user has access
    %to - Assume we have access to them all but we aren't sure of their
    %location, this will trigger the 'uigetfile' / 'uigetdir' calls.
    properties
        %Access to the FAME output folder (.pv8)
        HasPV8    = true;
        %Access to the FAME geometry file (.tux)
        HasTuxXML = true
        %Access to the FAME weights files (.205)
        Has205    = true;
        %Access to the FAME fuel distribution files
        HasFuel   = true;
    end
    
    % These properties are not seen on command line but can be accessed
    properties(Hidden = true)
        Info            % General information of number of nodes etc.
        Inp             % All the files that were parsed
        Fame            % FAME results
    end
    
    %% Methods
    methods
        
        % Contructor calling different classes
        function obj = Fame2mat(~)
            obj.Dirs = Dirs;
            obj.Opts = Options;
            obj.Mdl  = BeamModel;
            obj.Inp  = Inputfiles;
        end
        
        %% Run process to extract data and build NEOCASS and MSC.ADAMS models
        % Extract data from FAME model and construct NEOCASS and MSC.ADAMS
        % input files
        function obj = run(obj, logfcn)
            
            % N.B char(10) is a newline
            
            if nargin < 2
                
                logfcn = @(s) fprintf('%s\n', s);
                
                %Make some pretty display if we are printing to the command
                %window
                logfcn(' ===============================================================================');
                logfcn(char(10)); %#ok<*CHARTEN>
                logfcn('                             CONVERT FAME MODEL');
                logfcn(char(10))
                logfcn(sprintf('   USER NAME :   %s \n', upper(getenv('username'))))
                logfcn(sprintf('        DATE :   %s \n', datestr(now)))
                logfcn(char(10))
                logfcn(' -------------------------------------------------------------------------------');
                
            else
                
                logfcn('Converting FAME beam model');
                
            end
            
            try
                
                obj = getFameWing(obj, logfcn);
                obj = getFameResults(obj, logfcn);
                if ~isempty(obj.Inp.aeroPlanformFiles{1})
                    obj = getFameAero(obj, logfcn);
                    obj = getFameStructure(obj, logfcn);
                    obj = remeshStructure(obj, logfcn);
                    obj = defineEmpStructure(obj); %No fprintf statements
                    obj = setRbe0Nodes(obj, logfcn);
                    obj = setControlSurfaces(obj, logfcn);
                    obj = setAeroStructureSplines(obj, logfcn);
                    obj = getMassData(obj, logfcn);
                    obj = getFameEngine(obj); %No fprintf statements
                    obj = getFameAirfoils(obj, logfcn);
                    obj = setReflectWing(obj, logfcn);
                    obj = setInfo(obj); %No fprintf statements
                    obj = writeNeocassFiles(obj);
                    obj = writeSizingFile(obj);
                    obj = writeAdamsFiles(obj);
                end
                
                if nargin < 2
                    logfcn(char(10));
                    logfcn('Finished!');
                    logfcn('===============================================================================');
                else
                    logfcn('FAME beam model successfully converted');
                end
                
            catch runTimeError
                %                 save('ErrorLog','obj','runTimeError');
                rethrow(runTimeError);
            end
            
        end
        
        function obj = getFameWing(obj, logfcn)
            
            if nargin < 2
                logfcn = @(s) fprintf('%s\n', s);
            end
            
            logfcn('Reading FAME input file ... ');
            
            % Read the main .fm4 input file and extract as much data from
            % there as possible.
            
            if isempty(obj.Inp.fm4File{1}) 
                
                [fName,pName,idx] = uigetfile([pwd '*.fm4'],'Pick the fm4 file');
                if idx~=0
                    obj.Inp.fm4File{1} = [pName,fName];
                else
                    error('Unable to proceed without an *.fm4 file');
                end
            end
            
            obj.Dirs.fameFolder = fileparts(fileparts(obj.Inp.fm4File{1}));
            
            % Check the fame folder
            fs = dir(obj.Dirs.fameFolder);
            filecheck = sum(ismember({fs.name},{'__mirror','export','1_structure'}));
            
            if filecheck ~=3 && obj.HasPV8
                
                obj.Dirs.fameFolder = uigetdir(obj.Dirs.fameFolder,'Pick the FAME directory e.g. *-pv8');
                
            end
            
            FameFM4Structure = getFameInput(obj.Inp.fm4File{1});
            
            % Locate the tux_xml file:
            
            filelist = dir(fullfile([obj.Dirs.TUXFolder,filesep,'*.tux_xml']));
            
            if isempty(filelist) && obj.HasTuxXML
                [fName,pName,idx] = uigetfile([obj.Inp.fm4File{1} '*xml'],'Select the TUX_xml file');
                if idx~=0
                    obj.Inp.TUXFile{1} = [pName,fName];
                else
                    obj.Inp.TUXFile{1} = [];
                    logfcn('**WARNING** Proceeding without the TUX XML file. ');
                    %error('Unable to proceed without an *.tux_xml file');
                end
            else
                obj.Inp.TUXFile{1} = fullfile(filelist.folder,filelist.name);
            end
            
            if ~isempty(obj.Inp.TUXFile{1})
                
                % Read the tux_xml file and pass into the FameTUXStrucutre
                FameTUXStructure = xml2structFAME(obj.Inp.TUXFile{1});
                obj.Fame.TUXInput = FameTUXStructure;
            else
                FameTUXStructure = [];                
            end

            % Store fm4 and TUX structures somewhere in the object
            obj.Fame.fm4Input = FameFM4Structure;
            
            % Load some of that data into a more organised geometry
            % property
            obj.Fame.Geometry = sortFameInput(FameFM4Structure,logfcn);
            
            % Force the code to not deal with the empennage if it has found
            % no information about the tail
            
            if isempty(obj.Fame.fm4Input.HORIZONTAL_TAILPLANE.GLOBAL_GEOMETRY) && isempty(FameTUXStructure)
                obj.Opts.Geom.addTail = 0;
            end
            
            % Extract the WingOutline and store that in the geometry field
            folderName  = obj.Dirs.fameFolder;
            pathName    = ['export',filesep,'fmwpp01',filesep,'flexible'];
            
            % Check to see whether this folder exists
            if isdir([folderName,filesep,pathName])
                
                obj.Inp.aeroPlanformFiles{1} = fullfile(folderName,pathName,'wng_planform.fmwpp01');
                obj.Inp.aeroPlanformFiles{2} = fullfile(folderName,pathName,'htp_planform.fmwpp01');
                
                obj.Fame.Geometry.Wing.Outline = getWingOutline(obj.Inp.aeroPlanformFiles{1});
                
                if obj.Opts.Geom.addTail
                    
                    if ~isempty(FameTUXStructure)
                    
                    % Define the HTP and VTP outlines based on the tux_xml file
                    HTP_ref_point = cell2mat(textscan(FameTUXStructure.airplane.htp.AxisOffset.Attributes.value,'%f','Delimiter',','));
                    obj.Fame.Geometry.HTP.ref_point = HTP_ref_point';
                    HTP_y_le = cell2mat(textscan(FameTUXStructure.airplane.htp.planform.yle.Attributes.value,'%f','Delimiter',','));
                    HTP_y_te = cell2mat(textscan(FameTUXStructure.airplane.htp.planform.yte.Attributes.value,'%f','Delimiter',','));
                    HTP_x_le = cell2mat(textscan(FameTUXStructure.airplane.htp.planform.xle.Attributes.value,'%f','Delimiter',','));
                    HTP_x_te = cell2mat(textscan(FameTUXStructure.airplane.htp.planform.xte.Attributes.value,'%f','Delimiter',','));
                    HTP_z_te = zeros(size(HTP_x_te));
                    HTP_z_le = zeros(size(HTP_x_te));
                    
                    obj.Fame.Geometry.HTP.Outline = HTP_ref_point' + [HTP_x_le,HTP_y_le,HTP_z_le;flipud(HTP_x_te),flipud(HTP_y_te),flipud(HTP_z_te)];
                    
                    obj.Fame.Geometry.HTP.ControlSurfaces.elevator.eta_beg_le  = 0.05;
                    obj.Fame.Geometry.HTP.ControlSurfaces.elevator.eta_end_le  = 1;
                    obj.Fame.Geometry.HTP.ControlSurfaces.elevator.eta_beg_te  = 0.05;
                    obj.Fame.Geometry.HTP.ControlSurfaces.elevator.eta_end_te  = 1;
                    obj.Fame.Geometry.HTP.ControlSurfaces.elevator.chrd_beg_le = 0.65;
                    obj.Fame.Geometry.HTP.ControlSurfaces.elevator.chrd_end_le = 0.65;
                    
                    VTP_ref_point = cell2mat(textscan(FameTUXStructure.airplane.vtp.AxisOffset.Attributes.value,'%f','Delimiter',','));
                    VTP_y_le = cell2mat(textscan(FameTUXStructure.airplane.vtp.planform.yle.Attributes.value,'%f','Delimiter',','));
                    VTP_y_te = cell2mat(textscan(FameTUXStructure.airplane.vtp.planform.yte.Attributes.value,'%f','Delimiter',','));
                    VTP_x_le = cell2mat(textscan(FameTUXStructure.airplane.vtp.planform.xle.Attributes.value,'%f','Delimiter',','));
                    VTP_x_te = cell2mat(textscan(FameTUXStructure.airplane.vtp.planform.xte.Attributes.value,'%f','Delimiter',','));
                    VTP_z_te = zeros(size(VTP_x_te));
                    VTP_z_le = zeros(size(VTP_x_te));
                    
                    obj.Fame.Geometry.VTP.Outline = VTP_ref_point' + [VTP_x_le,VTP_z_le,VTP_y_le;flipud(VTP_x_te),flipud(VTP_z_te),flipud(VTP_y_te)];
                    obj.Fame.Geometry.VTP.ref_point = VTP_ref_point';
                    obj.Fame.Geometry.VTP.span      = max(VTP_y_te);
                    obj.Fame.Geometry.VTP.ControlSurfaces.rudder.eta_beg_le = 0.05;
                    obj.Fame.Geometry.VTP.ControlSurfaces.rudder.eta_end_le = 0.7;
                    obj.Fame.Geometry.VTP.ControlSurfaces.rudder.eta_beg_te = 0.05;
                    obj.Fame.Geometry.VTP.ControlSurfaces.rudder.eta_end_te = 0.7;
                    obj.Fame.Geometry.VTP.ControlSurfaces.rudder.chrd_beg_le = 0.65;
                    obj.Fame.Geometry.VTP.ControlSurfaces.rudder.chrd_end_le = 0.65;
                    
                    else
                        
                        folderName  = obj.Dirs.fameFolder;
                        pathName    = ['export',filesep,'fmwpp01',filesep,'flexible'];
                        fileName    = 'htp_planform.fmwpp01';
                        
                        HTPfile = fullfile(folderName,pathName,fileName);
                        
                        obj.Fame.Geometry.HTP.Outline = getHTPOutline(HTPfile);
                        
                        obj.Fame.Geometry.HTP.ControlSurfaces.elevator.eta_beg_le  = 0.05;
                        obj.Fame.Geometry.HTP.ControlSurfaces.elevator.eta_end_le  = 1;
                        obj.Fame.Geometry.HTP.ControlSurfaces.elevator.eta_beg_te  = 0.05;
                        obj.Fame.Geometry.HTP.ControlSurfaces.elevator.eta_end_te  = 1;
                        obj.Fame.Geometry.HTP.ControlSurfaces.elevator.chrd_beg_le = 0.65;
                        obj.Fame.Geometry.HTP.ControlSurfaces.elevator.chrd_end_le = 0.65;
                        
                    end
                else
                    if isfield(obj.Fame.Geometry,'HTP')
                        obj.Fame.Geometry = rmfield(obj.Fame.Geometry,'HTP');
                    end
                end
            else
                obj.Inp.aeroPlanformFiles{1} = [];
                obj.Inp.aeroPlanformFiles{2} = [];
                if isfield(obj.Fame.Geometry,'HTP')
                    obj.Fame.Geometry = rmfield(obj.Fame.Geometry,'HTP');
                end
            end
            
            logfcn('FAME input file reading complete.');
        end
        
        %% Obtain FAME data
        % Extract FAME results and place into model object
        function obj = getFameResults(obj, logfcn)
            
            if nargin < 2
                logfcn = @(s) fprintf('%s\n', s);
            end
            
            % Determine whether the folder in the Dirs object is ok.
            FolderFlag = 0;
            
            if isempty(obj.Dirs.fameFolder)
                warning('No FAME folder exists in the "Dirs" object');
                %                 error('FAME2MAT:DirNotExist','Fame directory not defined. Please add a directory in the "Dirs" object');
                FolderFlag = 1;
            elseif obj.Dirs.fameFolder == 0
                warning('The folder in the "Dirs" object does not exist');
                FolderFlag = 1;
            elseif ~exist(obj.Dirs.fameFolder)
                warning('The folder in the "Dirs" object does not exist');
                FolderFlag = 1;
            end
            
            if FolderFlag == 1
                
                FAME_folder = uigetdir([],'Choose FAME directory');
                if ischar(FAME_folder)
                set(obj.Dirs.famefolder,FAME_folder);
                end
                
            end
            
            % TODO Request the user for where the results files are
            
            logfcn('Reading FAME result files ... ');
            
            obj = getFameOutput(obj, logfcn);
            
            logfcn('FAME result file reading complete.');
            
            % Read VLM mesh
            
        end
        
        %% Extract FE data
        % Extract NASTRAN data from file generated by FAME. All data
        % placed in objects under the |obj.Mdl| object.
        function obj = getFameStructure(obj, logfcn)
            
            if nargin < 2
                logfcn = @(s) fprintf('%s\n', s);
            end
            
            if isempty(obj.Dirs.fameFolder)
                error('FAME2MAT:DirNotExist','Fame directory not defined. Please add a directory in the "Dirs" object');
            end
            
            logfcn('Reading FAME structural data files ...');
            obj = getNastran(obj, logfcn);
            logfcn('FAME structural data file reading complete.');
        end
        
        %% Remesh the beam model
        % The FAME FE models can sometimes contain a very fine mesh, hence
        % we use this method to remesh the structure. We check if the
        % |remeshStruct| structure of the |Struct| object is empty.
        function obj = remeshStructure(obj, logfcn)
            
            if nargin < 2
                logfcn = @(s) fprintf('%s\n', s);
            end
            
            if ~isempty(obj.Opts.Struct.remeshStruct.eta)
                logfcn('Remeshing structure ...');
                obj = remeshstruct2(obj);
                logfcn('Structure successfully remeshed.');
            end
            
        end
        
        function obj = defineEmpStructure(obj)
            
            if obj.Opts.Geom.addTail
                obj = generateEmpStructure(obj);
            end
        end
        
        %% Merge mass data
        % Merge mass data from FAME into |Conm2| objects. The .205 file is
        % parsed and then the mass data is divided into strips that are
        % either parallel to the XZ-plane or to the elastic axis of the
        % beam model. A polygon is drawn that incorporates the cut planes
        % and the leading and trailing edges. The points inside the polygon
        % are then determined, and the mass around the grid point inside
        % the polygin is calculated. The mass is written as a |Conm2| card
        % into the NASTRAN file.
        function obj = getMassData(obj, logfcn)
            
            if nargin < 2
                logfcn = @(s) fprintf('%s\n', s);
            end
            
            logfcn('Merging MASS data ...');
            
            % Merge structural and fuel data from FAME into Conm2 objects
            %sobj = mergeStructMass(obj, logfcn);
            % Check if the folder exists;
            findfuelpath = isdir([obj.Dirs.fameFolder,filesep,obj.Dirs.fameFuelFolder]);
            
            if findfuelpath == 1
                obj.Dirs.fameFuelFolder = [obj.Dirs.fameFolder,filesep,obj.Dirs.fameFuelFolder];
            end
            
            obj = mergeStructMass(obj, logfcn);
            obj = mergeFuelMass(obj, logfcn);
            
            obj = defineextramasses(obj);
            
            logfcn('MASS data successfully merged.');
            
        end
        
        %% Obtain FAME aerodynamic data
        % Read aerodynamic definitions data from FAME model
        function obj = getFameAero(obj, logfcn)
            
            if nargin < 2
                logfcn = @(s) fprintf('%s\n', s);
            end
            
            if isempty(obj.Dirs.fameFolder)
                error('FAME2MAT:DirNotExist','Fame directory not defined. Please add a directory in the "Dirs" object');
            elseif isempty(obj.Fame)
                error('FAME2MAT:ResultsEmpty','Fame results empty. Run "getFameStructure" method first');
            end
            
            logfcn('Reading FAME aerodynamic definition files ...');
            
            % Read and store wing geometry
            obj = fameAero2Obj(obj,'Wing');
            
            if obj.Opts.Geom.addTail
                
                if isfield(obj.Fame.Geometry,'HTP')
                    % Read and store HTP geometry
                    obj = fameAero2Obj(obj,'HTP');
                end
                
                if isfield(obj.Fame.Geometry,'VTP')
                    obj = fameAero2Obj(obj,'VTP');
                end
                
                % Read and store VTP geometry
                %obj = fameAeroConverterVTP(obj);
                
            end
            
            logfcn('FAME aerodynamic definition file reading complete.');
            
        end
        
        %% Aerodynamic nodes
        % Set the position of the aerodynamic nodes. In this case nodes are
        % created parallel to the XZ-plane or perpendicular to the elastic
        % axis, on the leading and trailing edges. These nodes are used for
        % the wing shape definition when interpolating the aero and
        % structural meshes.
        function obj = setRbe0Nodes(obj, logfcn)
            
            if nargin < 2
                logfcn = @(s) fprintf('%s\n', s);
            end
            
            logfcn('Creating RBE0 defintition ...');
            
            obj.Mdl.Rbe0 = Rbe0;
            
            % Create all the aero nodes
            obj = rbe0Generator(obj,'Wing');
            
            if obj.Opts.Geom.addTail
                if isfield(obj.Fame.Geometry,'HTP')
                    obj = rbe0Generator(obj,'HTP');
                    %obj = rbe0GeneratorHTP(obj);
                end
                if isfield(obj.Fame.Geometry,'VTP')
                    obj = rbe0Generator(obj,'VTP');
                    %obj = rbe0GeneratorVTP(obj);
                end
            end
            
            logfcn('RBE0 defintition complete.');
            
        end
        
        %% Controls surfaces
        % Define the control surface positions. This information is
        % contained in the |Fame| object. At the moment we are only
        % including ailerons.
        function obj = setControlSurfaces(obj, logfcn)
            
            if nargin < 2
                logfcn = @(s) fprintf('%s\n', s);
            end
            
            logfcn('Set control surfaces definition ...');
            
            if obj.Opts.Aero.wing_addControls == 1
                
                obj = controlSurfaceGenerator(obj,'Wing');
                
                if obj.Opts.Geom.addTail
                    if isfield(obj.Fame.Geometry,'HTP')
                        obj = controlSurfaceGenerator(obj,'HTP');
                    end
                    if isfield(obj.Fame.Geometry,'VTP')
                        obj = controlSurfaceGenerator(obj,'VTP');
                    end
                end
            end
            
            logfcn('Control surfaces defined.');
        end
        
        %% Connect aerodynamic and structural elements
        % Define the spline that connects the aerodynamic and structural
        % definitions
        function obj = setAeroStructureSplines(obj, logfcn)
            
            if nargin < 2
                logfcn = @(s) fprintf('%s\n', s);
            end
            
            logfcn('Defining aerodynamic and structural splines ...');
            
            obj.Mdl.Spline = Spline;
            obj.Mdl.Sets   = Sets;
            % Allocate the interpolation sets and spline cards
            obj = writeSpline(obj,'Wing');
            
            if obj.Opts.Geom.addTail
                
                if isfield(obj.Fame.Geometry,'HTP')
                    obj = writeSpline(obj,'HTP');
                    %obj = writeSplineHTP(obj,'HTP');
                end
                
                if isfield(obj.Fame.Geometry,'VTP')
                    obj = writeSpline(obj,'VTP');
                    %obj = writeSplineVTP(obj,'HTP');
                end
            end
            
            logfcn('Aerodynamic and structural splines defined.');
            
        end
        
        %% Engine defintion
        % Obtain engine mass and location and connect with closest node via
        % |Rbe2| element.
        function obj = getFameEngine(obj)
            
            if obj.Opts.Struct.addEngine==1
                %fprintf(' Reading FAME engine data files ...\n');
                
                %obj = getFameEngineDataFcn(obj);
            end
            
        end
        
        %% Airfoil Definitions
        % Obtain engine mass and location and connect with closest node via
        % |Rbe2| element.
        function obj = getFameAirfoils(obj, logfcn)
            
            if nargin < 2
                logfcn = @(s) fprintf('%s\n', s);
            end
            
            
            if obj.Opts.Aero.wing_addCamber
                logfcn('Creating Airfoil data files ...');
                obj = interpprofneocassformat(obj);
                logfcn('Airfoil data files successfully created.');
            end
            
        end
        
        %% Reflect the wing
        % Create mirror image of wing structure and aerodynamic definitions
        function obj = setReflectWing(obj, logfcn)
            
            if nargin < 2
                logfcn = @(s) fprintf('%s\n', s);
            end
            
            if obj.Opts.Geom.reflect==1
                
                logfcn('Reflecting wing definition ...');
                
                % Do structure mirroring first
                obj = symWingStruct(obj);
                
                % Now do aero mirroring
                obj = symWingAero(obj);
                
                % Reflect thrust
                obj = symWingThrust(obj);
                
                % Reflect SPC's
                obj = symWingSpc(obj);
                
                logfcn('Wing definition successfully reflected.');
            end
            
        end
        
        %% Information about model
        % Set accounting of number of different elements in information
        % object. Contains information about the number of elements of each
        % object.
        function obj = setInfo(obj)
            
            obj.Mdl.Info.nGrid    = numel(obj.Mdl.Grid);
            obj.Mdl.Info.nMgrid   = [];
            obj.Mdl.Info.nPbar    = numel(obj.Mdl.Pbar);
            obj.Mdl.Info.nBar     = numel(obj.Mdl.Cbar);
            obj.Mdl.Info.nMat     = numel(obj.Mdl.Mat);
            obj.Mdl.Info.nForce   = numel(obj.Mdl.Force);
            obj.Mdl.Info.nMoment  = numel(obj.Mdl.Moment);
            obj.Mdl.Info.nConm2   = numel(obj.Mdl.Conm2);
            obj.Mdl.Info.nCaero   = numel(obj.Mdl.Caero);
            obj.Mdl.Info.nRbe0    = numel(obj.Mdl.Rbe0);
            obj.Mdl.Info.nSgrid   = [];
            obj.Mdl.Info.nSets    = numel(obj.Mdl.Sets);
            obj.Mdl.Info.nEngines = numel(find([obj.Mdl.Grid.id]>300000));
            obj.Mdl.Info.nRbe2    = numel(obj.Mdl.Rbe2);
            obj.Mdl.Info.nThrust  = numel(obj.Mdl.Thrust);
            obj.Mdl.Info.nSpc     = numel(obj.Mdl.Spc);
            obj.Mdl.Info.nLoads   = [];
            
        end
        
        %% Plot NEOCASS model
        % Plot stick model of the NEOCASS structural model
        function obj = plotNeocassModel(obj)
            
            plotNeocassModelFunc(obj);
            
        end
        
        %% Create NEOCASS files
        % Write the output files that are neeed by NEOCASS. Format similar
        % to NASTRAN files.
        function obj = writeNeocassFiles(obj)
            
            if obj.Opts.Outputs.writeNeocassFiles == false
                return
            end
            
            if ~exist(obj.Dirs.writeFolder,'dir')
                mkdir(obj.Dirs.writeFolder)
            end
            
            fprintf(' Writing NEOCASS input files ...\n');
            
            % Write data to a main structure-aerodynamic input file
            
            if ~exist([obj.Dirs.writeFolder filesep 'Neo_model'],'dir')
                mkdir([obj.Dirs.writeFolder filesep 'Neo_model'])
            end
            fame2neo(fullfile(obj.Dirs.writeFolder,'Neo_model','FameModel_Neo.dat'),obj);
            
            % Generate interpolated airfoil profiles
            starty=[];
            for i = 1:numel(obj.Mdl.Caero)
                starty = [starty,obj.Mdl.Caero(1,i).startY];
            end
            
            endy=[];
            for i = 1:numel(obj.Mdl.Caero)
                endy = [endy,obj.Mdl.Caero(1,i).startY + obj.Mdl.Caero(1,i).b];
            end
            
            y_loc = [starty,endy(end)];
            y_loc = y_loc/y_loc(end);
            
            % Generate fuel decks
            fuel2neo(obj.Dirs.writeFolder,obj);
            
            % Generate Payload files
            writePayloadFiles(obj.Dirs.writeFolder,obj);
            
            % Write TRIM input file
            writeTrimInputfile(obj.Dirs.writeFolder,obj);
            
            % Write MAIN INPUT FILES
            writeMainInputfile(obj.Dirs.writeFolder,obj);
            
            % Write an input file with the correct thicknesses
            writeThicknessInput(obj.Dirs.writeFolder,obj);
        end
        
        %% Create NEOCASS files
        % Write the output files that are neeed by NEOCASS. Format similar
        % to NASTRAN files.
        function obj = writeSizingFile(obj)
            
            for idx = 1:numel(obj.Fame.LoadCases)
                obj.Fame.LoadCases(idx).PayloadMass = obj.Fame.LoadCases(idx).AC_mass - 2*obj.Fame.LoadCases(idx).fuel - obj.Fame.Geometry.Weights.oew;
                obj.Fame.LoadCases(idx).PayloadMass = [];

            end
            
            if obj.Opts.Outputs.writeNeocassFiles == false
                return
            end
            obj = fame2LoadCase(obj.Dirs.writeFolder,obj);
            
        end
        
        %% Create MSC.ADAMS file
        % A beam model is created from the NASTRAN data, with accompanying
        % mass data, as well as panels for the aerodynamic panels.
        function obj = writeAdamsFiles(obj)
            % Write the output files for ADAMS
            if obj.Opts.Outputs.writeAdamsFiles == false
                return
            end
            
            if ~exist(obj.Dirs.writeFolder,'dir')
                mkdir(obj.Dirs.writeFolder)
            end
            
            fprintf(' Writing MSC.ADAMS command file ...\n');
            
            % Write data to a main structure-aerodynamic input file
            fame2adams(fullfile(obj.Dirs.writeFolder,'Fame2AdamsConversion.cmd'),obj);
            
        end
        
        %% Set field values of |model| object
        % Define a set method for this class to make sure relevant
        % fields are set and whole object not redefined in a similar
        % way to structure behaviour
        function obj = set.Mdl(obj,Value)
            
            if isa(Value,'BeamModel')
                obj.Mdl = Value;
            elseif isobject(Value) || isstruct(Value)
                fieldsStruct = fieldnames(Value);
                fieldsObject = fieldnames(obj.Mdl);
                fieldsObjectReduced = regexprep(fieldsObject,'Obj','');
                
                [~,iVal,iObj] = intersect(fieldsStruct,fieldsObjectReduced);
                
                for i = 1:numel(iObj)
                    obj.Mdl.(fieldsObject{iObj(i)}) = Value.(fieldsStruct{iVal(i)});
                end
            end
            
        end
    end
end
