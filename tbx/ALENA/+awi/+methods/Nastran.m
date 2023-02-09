classdef Nastran < awi.methods.Analysis
    %Nastran Generic analysis class containing methods and properties
    %related to MSC.Nastran analysis.
    
    properties (Access = public)
        %Path to the MSC.Nastran executable
        NastranExe
        %Reference node for applying SPC and SUPORT entrie
        RefNode
        WingNode
        FWTNode
    end
    
    properties (Constant, Hidden = true)
        %Maximum amount of time (in seconds) to wait for the .f06 file
        F06MaxWait = 60;
        %Invalid characters in filenames - Replacd by '_'
        InvalidFileChar = {'.', '@', '(', ')', '{', '}', '[', ']', '__'};
        %Initial ID number for any extra bulk data not provided in the
        %AnalysisModel.
        BulkDataID0 = 250;
    end
    
    methods % set / get
        function set.NastranExe(obj, val)  %set.NastranExe
            %set.NastranExe Set method for the property 'NastranExe'.
            validateattributes(val, {'char'}, {'row', 'nonempty'}, ...
                class(obj), 'NastranExe');
            obj.NastranExe = val;
        end
        function val = get.NastranExe(obj) %get.NastranExe
            %get.NastranExe Get method for the property 'NastranExe'.
            %
            % If 'NastranExe' has not been initialized then use 'getpref'.
            
            val = obj.NastranExe;
            if isempty(val)
                %Use 'setpref'/'getpref'
                if ispref('ALENA_Nastran', 'nastran_exe')
                    val = getpref('ALENA_Nastran', 'nastran_exe');
                else
                    %Ask the user
                    [name, path] = uigetfile({'*.exe', 'Executable File (*.exe)'}, ...
                        ['Select the MSC.Nastran executable file ', ...
                        '(e.g. \...\nastranw.exe)']);
                    %Check the output
                    if isnumeric(name) || isnumeric(path)
                        warndlg(['Unable to run analysis as the '  , ...
                            'path to the MSC.Nastran executable '  , ...
                            'has not been set. Update your user '  , ...
                            'preferences and re-run the analysis.'], ...
                            'Warning - Unable to run analysis', 'modal');
                        return
                    end
                    %Check the path for spaces - Enclose with ""
                    val = fullfile(path, name);
                    if contains(val,' ')
                        val = ['"',val,'"'];
                    end
                    %Update the preferences
                    setpref('ALENA_Nastran', 'nastran_exe', val);
                end
            end
        end
        function set.RefNode(obj, val)
            %set.RefNode Set method for the property 'RefNode'.
            
            validateattributes(val, {'awi.fe.Node'}, {'scalar'}, ...
                class(obj), 'RefNode');
            obj.RefNode = val;
            
        end
        function set.WingNode(obj, val)
            %set.WingNode Set method for the property 'WingNode'.
            
            validateattributes(val, {'awi.fe.Node'}, {'scalar'}, ...
                class(obj), 'WingNode');
            obj.WingNode = val;
            
        end
        function set.FWTNode(obj, val)
            %set.FWTNode Set method for the property 'FWTNode'.
            
            validateattributes(val, {'awi.fe.Node'}, {'scalar'}, ...
                class(obj), 'FWTNode');
            obj.FWTNode = val;
            
        end
        function val = get.RefNode(obj)
            %get.RefNode Get method for the property 'RefNode'.
            
            val = obj.RefNode;
            if isempty(val) && ~isempty(obj.AnalysisModel)
                %Where will the constraints be applied?
                %   - Just go with the start of the fuselage for now
                AllFEM = flatlist(obj.AnalysisModel);
                Fuselage = AllFEM(ismember({AllFEM.Name}, 'Fuselage'));
                if isempty(Fuselage)
                    for i = 1:length(AllFEM)
                        if ~isempty(AllFEM(i).Nodes)
                            val = AllFEM(i).Nodes(1);
                            break;
                        end
                    end
                else
                    num=numel(Fuselage.Nodes)/5;
                    val  = Fuselage.Nodes(1:num);
                    %select all along the fuselage (charles)
                end
            end
            
        end
        % select wing nodes on FE model(charles)
        function val = get.WingNode(obj)
            %get.RefNode Get method for the property 'WingNode'.
            
            val = obj.WingNode;
            if isempty(val) && ~isempty(obj.AnalysisModel)
                AllFEM = flatlist(obj.AnalysisModel);
                Wingbox_right = AllFEM(ismember({AllFEM.Name}, 'A320Wing_right'));
                if ~isempty(Wingbox_right)
                    num=(numel(Wingbox_right.Nodes))/3;
                    val  = Wingbox_right.Nodes(1:num);
                end
            end
            
        end
        
        % select all nodes on the FWT(charles)
        function val = get.FWTNode(obj)
            val = obj.FWTNode;
            if isempty(val) && ~isempty(obj.AnalysisModel)
                AllFEM = flatlist(obj.AnalysisModel);
                FWT = AllFEM(ismember({AllFEM.Name}, 'A320Wing_right_FWT'));
                if ~isempty(FWT)
                    num=(numel(FWT.Nodes))/3;
                    val  = FWT.Nodes(1:num);
                end
            end
            
        end
    end
    
    methods % static analysis
        function StaticResults = static(obj, StaticOptions, resDir)
            %static Executes a SOL 101 or SOL 400 analysis in MSC.Nastran
            %using the analysis model 'obj.AnalysisModel' and the options
            %in 'StaticOptions'.
            
            StaticResults = [];
            
            if nargin < 3
                resDir = obj.makeDefaultAnalysisDirectory(analysisType);
            end
            assert(isfolder(resDir), ['Expected the output directory ', ...
                'to be a valid folder location.']);
            
            %Pass it on
            static@awi.methods.Analysis(obj, StaticOptions);
            
            %Write the header file
            switch StaticOptions.StructuralResponse
                case 'linear'
                    staticFile = obj.writeLinearStaticFile(StaticOptions, resDir);
                case 'nonlinear'
                    staticFile = obj.writeNonlinearStaticFile(StaticOptions, resDir);
            end
            %staticFile = fullfile(resDir, 'nonlinear_static_analysis.dat');
            
            %Run the analysis
            %   - N.B. Will change the current directory to force Nastran
            %     to output the files in the correct location
%             obj.runNastran(staticFile, 'MonitorProgress', 'Nonlinear');
            obj.runNastran(staticFile, 'MonitorProgress', 'Default');
            
            %Extract the results from the Nastran output
            NastranResults = obj.extractNastranResults(staticFile);
            
            %Generate the 'awi' results objects
            StaticResults = obj.generateAWIResults(NastranResults);
            
        end
        function staticFile = writeStaticFile(obj, StaticOptions, outDir)
            
            %Update the flags of any loads to 'global' instead of local
            
        end
        function runFile = writeNonlinearStaticFile(obj, StaticOptions, outDir, varargin)
            %writeNonlinearStaticFile Writes the header file for a SOL 400
            %nonlinear static analysis.
            
            %Grab the FEM
            FEM = obj.AnalysisModel;
            
            %Parse
            p = inputParser;
            addParameter(p, 'DatFilename', [], @(x)validateattributes(x, {'char'}, {'row'}));
            parse(p, varargin{:});
            if nargin < 3
                outDir = obj.makeDefaultAnalysisDirectory(analysisType);
            end
            assert(isfolder(outDir), ['Expected the output directory ', ...
                'to be a valid folder location.']);
            
            %Define additional data for nonlinear analysis
            NonlinFEM = awi.fe.FEModel;
            NonlinFEM.Name = 'NonlinearModelData';
            
            %Define additional nodes for specifiying orientation of
            %follower loads
            bFollower = i_defineFollowerNodes(FEM, NonlinFEM);
            
            function bFollower = i_defineFollowerNodes(FEM, NonlinFEM)
                %i_defineFollowerNodes Defines a series of RBE spider nodes
                %for specifying the orientation of the follower loads and
                %assigns this information to the 'awi.fe.FEModel' object
                %'NonlinFEM'.
                
                bFollower = false;
                
                Loads = FEM.PointLoads;
                
                %Do we need to do anything?
                if isempty(Loads)
                    return
                end
                idx = ismember({Loads.LoadBehaviour}, 'follower');
                if ~any(idx)
                    return
                end
                bFollower = true;
            
                %Make a new node and RBE for each load
                Loads = Loads(idx);
                nLoad = numel(Loads);                
                OrientNodes = arrayfun(@(~) awi.fe.Node    , 1 : nLoad);
                SlaveRBE    = arrayfun(@(~) awi.fe.RigidBar, 1 : nLoad);
                
                %Grab load data
                LoadNodes = [Loads.Node];
                v = [Loads.Orientation];
                v = v ./ vecnorm(v);
                
                %Calculate position of nodes for defining follower loads
                idxL  = ismember({Loads.LoadCoordSys}, 'local');
                if any(idxL)
                    ind = find(idxL);
                    %   - Transform orientation of load into global CS
                    rotMatrix   = getRotationMatrix([LoadNodes(idxL).OutputCoordSys]);
                    v_ = arrayfun(@(i) rotMatrix(:, :, ind(i)) * v(:, i), 1 : nnz(idxL), 'Unif', false);
                    v(:, idxL) = horzcat(v_{:});
                end
                nodeX = [LoadNodes.X] + v;
                
                %Assign ID numbers to FE objects
                FEObj = [OrientNodes, SlaveRBE];
                FEM.assignIDnumbers(FEObj)
                id0 = FEM.MaxIDNumber;
                if isempty(id0)
                    id0 = max([FEObj.ID0]);
                end
                id = (id0 + 1 : id0 + numel(FEObj))';
                set(FEObj, {'ID'}, num2cell(id));
                
                
                
                %Assign data to objects & FEM
                set(Loads      , {'OrientationNode'}, num2cell(OrientNodes)');
                set(OrientNodes, {'X'}, num2cell(nodeX, 1)');
                set(SlaveRBE   , {'NodesI'}, num2cell(LoadNodes)');
                set(SlaveRBE   , {'NodesD'}, num2cell(OrientNodes)');
                set(SlaveRBE   , 'CN'      , '123456');
                addFEData(NonlinFEM, FEObj);
                
            end
            
            %Define the combination of all loads
            LoadSet = awi.fe.LoadSet;
            LoadSet.ID    = 11;
            LoadSet.Loads = FEM.PointLoads;
            LoadSet.Si    = ones(1, numel(LoadSet.Loads));
            
            %Export the FEM & the additional nonlinear information
            [~, includeFiles] = export(FEM, outDir, 'WriteHeaderFile', false);
            [~, nonlinFiles ] = export(NonlinFEM, outDir, 'WriteHeaderFile', false);
            includeFiles      = [includeFiles, nonlinFiles];
            
            %% Write the header file
            runFile = p.Results.DatFilename;
            if isempty(runFile)
                runFile = 'nonlinear_static_analysis';
            else
                runFile = obj.parseFilename(runFile);                
            end
            runFile = fullfile(outDir, [runFile, '.dat']);
            
            fid = fopen(runFile, 'w');
            
            %Executive Control
            awi.fe.FEBaseClass.writeHeading(fid, 'E X E C U T I V E  C O N T R O L');
            fprintf(fid, 'SOL %i\r\n', 400);
            fprintf(fid, 'ECHOOFF         $ SUPPRESSES THE ECHO OF EXECUTIVE CONTROL\r\n');
            fprintf(fid, 'CEND\r\n');
            
            %Case Control
            awi.fe.FEBaseClass.writeHeading(fid, 'C A S E  C O N T R O L');
            awi.fe.FEBaseClass.writeSubHeading(fid, 'O U T P U T  O P T I O N S');
            fprintf(fid, 'LINE = 99999   $ SPECIFIES THE NUMBER OF LINES PER PRINTED PAGE\r\n');
            fprintf(fid, 'ECHO = NONE    $ SUPPRESSES THE ECHO OF BULK DATA\r\n');
            %
            awi.fe.FEBaseClass.writeSubHeading(fid, 'O U T P U T  Q U A N T I T I E S');
            fprintf(fid, 'DISP  = ALL $ OUTPUT ALL DISPLACEMENTS\r\n');
            fprintf(fid, 'STRESS(VONMISES)  = ALL $ OUTPUT ALL STRESSES\r\n');
            fprintf(fid, 'FORCE  = ALL $ OUTPUT ALL FORCES AND MOMENTS\r\n');
            %
            awi.fe.FEBaseClass.writeSubHeading(fid, 'G L O B A L  C A R D S')
            fprintf(fid, 'ANALYSIS = NLSTATIC\n'); 
            fprintf(fid, 'SPC  = %i\r\n', 999);            
                        
            %Write STEP entries for nonlinear analysis
            fprintf(fid, 'STEP 1\n');
            fprintf(fid, '%-10s = %i\n', 'NLSTEP', 10);
            fprintf(fid, '%-10s = %i\n', 'LOAD'  , LoadSet.ID);
            
            %Bulk Data
            obj.defaultBulkStatement(fid);
            if bFollower
                fprintf(fid, 'PARAM,LGDISP,1\n');
            else
                fprintf(fid, 'PARAM,LGDISP,2\n');
            end
            %Write the cards
            awi.fe.FEBaseClass.writeColumnDelimiter(fid, 'normal');

            %   - SPC

            if isempty(FEM.Children.Children)

                fprintf(fid, '%-8s%-8i%-8i%-8i\r\n', 'SPC1', 999, 123456, FEM.Nodes(1).ID);

            else

                fprintf(fid, '%-8s%-8i%-8i%-8i%-8i%-8i%-8i%-8i%-8i%-8i%-8i\r\n', 'SPC1', 999, 123456, [FEM.Children(1).Children(1).Children.Nodes(1).ID, ...
                                                                                                       FEM.Children(1).Children(2).Children.Nodes(1).ID,...
                                                                                                       FEM.Children(1).Nodes(1).ID,...
                                                                                                       FEM.Children(1).Nodes(40).ID,...
                                                                                                       FEM.Children(1).Nodes(42).ID,...
                                                                                                       FEM.Children(1).Nodes(end/5).ID,...
                                                                                                       FEM.Children(1).Children(1).Nodes(4).ID,...
                                                                                                       FEM.Children(1).Children(2).Nodes(4).ID]);

            end

            %   - NLSTEP
            CTRLDEF = StaticOptions.CTRLDEF; % Option for controlling load iteration in SOL 400, {'', 'QLINEAR', 'MILDLY', 'SEVERELY'}
            IDAMP = 6; % Flag for activating artificial damping for static analysis for robust and stable analysis
            fprintf(fid, '%-8s%-8i%#-8.3f%-8s\n%-8s%-8s\n%-8s%-8s%-8i\n', ...
                'NLSTEP', 10, 1, CTRLDEF, blanks(8), 'ADAPT', blanks(8), blanks(8), IDAMP);
            
            % Adding gravity if stated in StaticOptions
            if StaticOptions.bModelGravity
                gVec = StaticOptions.GravityVector;
                gravID = LoadSet.LoadID(end)+1; 
                fprintf(fid, '%-8s%-8i%-8i%-8.5f%-8.3f%-8.3f%-8.3f\r\n', ...
                        'GRAV', gravID, 0, 1.0, gVec(1), gVec(2), gVec(3));
                % Add gravity ID to load set
                grav = awi.fe.PointLoad; % ... adding gravity as a "point load" because LoadSet only accepts point loads, create a new "awi.fe.Gravity" class?
                grav.SID = gravID;
              	LoadSet.Si = [LoadSet.Si 1]; % Adding gravity scale factor of 1
                LoadSet.Loads = [LoadSet.Loads grav]; % Add gravity to loads (unclear how LoadSet.LoadID gets written...)
            end
            
            %   - LOAD
            writeToFile(LoadSet, fid);
            
            %Additional bulk data files
            if ~iscell(includeFiles)
                includeFiles = {includeFiles};
            end
            awi.fe.FEBaseClass.writeSubHeading(fid, 'I N C L U D E  F I L E S');
            awi.fe.FEBaseClass.writeIncludeStatement(fid, includeFiles);
            
            %End of file
            fprintf(fid, 'ENDDATA\r\n');
            
            %Close the file
            fclose(fid);
            
        end
    end
    
    methods % trim analysis
        function TrimResults = trim(obj, Aircraft, LoadCases)
            %trim Executes a Solution 144 analysis in MSC.Nastran using the
            %analysis model 'obj.AnalysisModel', the aircraft data in
            %'Aircraft' and the load case information in 'LoadCases'.
            
            TrimResults = [];
            
            %Pass it on
            trim@awi.methods.Analysis(obj, Aircraft, LoadCases);
            
            %Write the header file
            trimFile = obj.writeTrimFile(Aircraft, LoadCases);
            
            %Run the analysis
            %   - Need to change the current directory to force Nastran to
            %     output the files in the correct location
            obj.runNastran(trimFile);
            
            %Extract the results from the Nastran output
            %   - The aero forces & pressures are in the .f06
            NastranResults = obj.extractNastranResults(trimFile, 'ReadF06', true);
            
            %Generate the 'awi.model.TrimResult' object
            TrimResults = obj.generateAWIResults(NastranResults);
            
        end
        function trimFile = writeTrimFile(obj, Aircraft, LoadCases, MassCases, outDir, varargin)
            %writeTrimFile Writes a header file for a static aeroelastic
            %analysis using SOL 144 for the flight point defined in the
            %'awi.model.LoadCase' object(s) 'LoadCases'.
            %
            % Multiple load cases will be defined over multiple subcases.
            
            p = inputParser;
            addParameter(p, 'DatFilename', [], @(x)validateattributes(x, {'char'}, {'row'}));
            p.addParameter('AoA', NaN);
            p.addParameter('ExtraCards',[]);
            parse(p, varargin{:});
            
            %Grab the FEM
            FEM = obj.AnalysisModel;
            
            %Check we can actually do a trim analysis
            obj.checkAircraftTrimData(Aircraft);
            
            %Grab the trim data
            ID0      = obj.BulkDataID0;
            TrimData = awi.methods.Nastran.getTrimData(FEM, Aircraft, LoadCases, MassCases, ID0, obj.RefNode,'AoA',p.Results.AoA);
            
            %% Set up the directory and print model bulk data
            
            %Create a folder for the output
            if nargin < 4
                dirName = [datestr(now, 30), '_trim_analysis'];
                [bSuccess, ~, ~] = mkdir(pwd, dirName);
                assert(bSuccess == true, ['Unable to create output ', ...
                    'directory for the trim analsysis.']);
                outDir = fullfile(pwd, dirName);
            end
            
            %Export the FEM
            [~, includeFiles] = export(FEM, outDir, 'WriteHeaderFile', false);
            
            %% Write the header file
            trimFile = p.Results.DatFilename;
            if ~isempty(trimFile)
                trimFile =  obj.parseFilename(trimFile);
            elseif numel(LoadCases) > 1 || isempty(LoadCases.Name)
                trimFile = 'trim_analysis';
            else
                trimFile = obj.parseFilename(LoadCases.Name);
            end
            trimFile = fullfile(outDir, [trimFile, '.dat']);
            
            fid = fopen(trimFile, 'w');
            
            %Executive Control
            awi.fe.FEBaseClass.writeHeading(fid, 'E X E C U T I V E  C O N T R O L');
            fprintf(fid, 'NASTRAN NLINES=999999 \r\n');
            
            % output Ajj and FFAJ matrix
            
            Ajj_file=strcat('''','AJJ.op4','''');
            fprintf(fid, '%s%s%s\r\n','ASSIGN output4=', Ajj_file ,',formatted,UNIT=11');
            
            FFAJ_file=strcat('''','FFAJ.op4','''');
            fprintf(fid, '%s%s%s\r\n','ASSIGN output4=', FFAJ_file ,',formatted,UNIT=12');
            
            fprintf(fid, 'SOL %i\r\n', 144);
            
            fprintf(fid, 'COMPILE PFAERO $\r\n');
            fprintf(fid, 'ALTER 278$ \r\n');
            fprintf(fid, 'OUTPUT4 AJJ,,,,//0/11///8 $ \r\n');
            fprintf(fid, 'COMPILE AESTATRS $ \r\n');
            
            ASDR=strcat('''','ASDR','''');
            fprintf(fid, '%s%s%s\r\n', 'ALTER', ASDR, '$' );
            fprintf(fid, 'OUTPUT4 FFAJ,,,,//0/12///8 $ \r\n');
            
            
            fprintf(fid, 'ECHOOFF         $ SUPPRESSES THE ECHO OF EXECUTIVE CONTROL\r\n');
            fprintf(fid, 'CEND\r\n');
            
            %Case Control
            awi.fe.FEBaseClass.writeHeading(fid, 'C A S E  C O N T R O L');
            awi.fe.FEBaseClass.writeSubHeading(fid, 'O U T P U T  O P T I O N S');
            fprintf(fid, 'LINE = 99999   $ SPECIFIES THE NUMBER OF LINES PER PRINTED PAGE\r\n');
            fprintf(fid, 'ECHO = NONE    $ SUPPRESSES THE ECHO OF BULK DATA\r\n');
            
            %symmetry conditions
            fprintf(fid, 'AECONFIG=AeroSG2D\r\n');
            fprintf(fid, 'AESYMXZ=Asymmetric\r\n');
%             fprintf(fid, 'AESYMXZ=Symmetric\r\n');
            fprintf(fid, 'AESYMXY=Asymmetric\r\n');
            
            %
            awi.fe.FEBaseClass.writeSubHeading(fid, 'O U T P U T  Q U A N T I T I E S');
            fprintf(fid, 'DISP(PLOT)  = ALL $ OUTPUT ALL DISPLACEMENTS\r\n');
            fprintf(fid, 'DISP = ALL $ OUTPUT ALL DISPLACEMENTS\r\n');
            fprintf(fid, 'FORCE(PLOT) = ALL $ OUTPUT ALL ELEMENT FORCES\r\n');
            fprintf(fid, 'AEROF       = ALL $ OUTPUT ALL AERO FORCES\r\n');
            fprintf(fid, 'APRESSURE   = ALL $ OUTPUT ALL AERO PRESSURES\r\n');
            fprintf(fid, 'STRESS   = ALL $ OUTPUT ALL STRESSES\r\n');
            fprintf(fid, 'FORCE   = ALL $ OUTPUT ALL FORCES AND MOMENT\r\n');
            fprintf(fid, 'TRIMF(APPLIED,AIR)   = ALL $ OUTPUT ALL TRIM LOAD\r\n');
            %
            awi.fe.FEBaseClass.writeSubHeading(fid, 'G L O B A L  C A R D S')
            fprintf(fid, 'SPC  = %i\r\n', TrimData.SPC_id);

            awi.fe.FEBaseClass.writeSubHeading(fid, 'S U B C A S E S')
            pad = blanks(2);
            for iS = 1 : numel(TrimData.Trim)
                fprintf(fid, 'SUBCASE %i\r\n', iS);
                fprintf(fid, '%sTRIM = %i\r\n', pad, TrimData.Trim(iS).SID);

                %Add GRAV card if the entry is not empty.
                if isfield(TrimData,'GRAV')

                    TrimData.LOAD_ID=TrimData.Trim(iS).SID*10;
                    fprintf(fid, '%sLOAD = %i\r\n', pad, TrimData.LOAD_ID);

                end
            end
            

            
            %Bulk Data
            obj.defaultBulkStatement(fid);
            %Write the cards
            mni.printing.bdf.writeColumnDelimiter(fid, 'normal');

            %   - GRAV
            if isfield(TrimData,'GRAV')
                mni.printing.cards.GRAV(TrimData.LOAD_ID,TrimData.GRAV,[0,0,-1],'CID',0).writeToFile(fid);
            end

            %   - SPC
            mni.printing.cards.SPC1(TrimData.SPC_id,...
                TrimData.SPCdof, [TrimData.RefGrid(46).ID]).writeToFile(fid);
            %   - SUPORT
            if ~isnan(TrimData.SUPdof)
                mni.printing.cards.SUPORT(TrimData.RefGrid(46).ID,...
                    TrimData.SUPdof).writeToFile(fid);
            end
            %   - AEROS
            obj.writeSteadyAeroEntry(fid, Aircraft);
%             obj.writeSteadyAeroEntry(fid, Aircraft,'SymXZ',true);

            
            %   - AESTAT
            aestat_fun = @(idx)mni.printing.cards.AESTAT(...
                TrimData.AESTAT_id(idx),TrimData.AESTAT{idx}).writeToFile(fid);
            arrayfun(aestat_fun,1:length(TrimData.AESTAT))
            %   - TRIM
            if LoadCases.Flexibility == 0
                mni.printing.cards.TRIM(TrimData.Trim.SID,...
                    TrimData.Trim.MACH,TrimData.Trim.Q,...
                    TrimData.Trim.TrimVar(:),...
                    'AEQR',0).writeToFile(fid);
            else
                mni.printing.cards.TRIM(TrimData.Trim.SID,...
                    TrimData.Trim.MACH,TrimData.Trim.Q,...
                    TrimData.Trim.TrimVar(:)).writeToFile(fid);
            end
            mni.printing.cards.PARAM('AUNITS','r',0.1020).writeToFile(fid)            
            % - Extra Cards
            arrayfun(@(x)x.writeToFile(fid),p.Results.ExtraCards)
            
            %Additional bulk data files
            if ~iscell(includeFiles)
                includeFiles = {includeFiles};
            end
            
            mni.printing.bdf.writeSubHeading(fid, 'I N C L U D E  F I L E S');
            includeFiles = cellfun(@(x)mni.printing.cards.INCLUDE(x),includeFiles);
            arrayfun(@(x)x.writeToFile(fid),includeFiles)
            
            %End of file
            fprintf(fid, 'ENDDATA\r\n');
            
            %Close the file
            fclose(fid);
            
        end
    end
    
    methods % gust analysis
        function GustResults = gust(obj, Aircraft, LoadCases)
            %gust Executes a Solution 146 analysis in MSC.Nastran using the
            %analysis model 'obj.AnalysisModel', the aircraft data in
            %'Aircraft' and the load case information in 'LoadCases'.
            
            GustResults = [];
            
%             for iG = 1 : numel(LoadCases)
            %Write the header file
            gustFile = obj.writeGustFile(Aircraft, LoadCases);
%         end

            %Run the analysis
            %   - Need to change the current directory to force Nastran to
            %     output the files in the correct location
            obj.runNastran(gustFile);
            
            %Extract the results from the Nastran output
            %   - The aero forces & pressures are in the .f06
%             NastranResults = obj.extractNastranResults(gustFile, 'ReadF06', true);
            
            %Generate the 'awi.model.TrimResult' object
%             TrimResults = obj.generateAWIResults(NastranResults);
            
        end
        function gustFile = writeGustFile(obj, Aircraft, LoadCases, MassCases, FlightPoints, outDir, varargin) %Added obj at the front
            %writeGustFile Writes a header file for a dynamic aeroelastic
            %analysis using SOL 146 for the flight point defined in the
            %'awi.model.LoadCase' object(s) 'LoadCases'.
            %
            % Multiple load cases will be defined over multiple subcases.
            
            if nargin < 3
                MassCases = [];
            end
            
            p = inputParser;
            addParameter(p, 'DatFilename', [], @(x)validateattributes(x, {'char'}, {'row'}));
            parse(p, varargin{:});
            
            %Grab the FEM
            FEM = obj.AnalysisModel;
            
            %Check we can actually do a trim analysis
            obj.checkAircraftTrimData(Aircraft);
            
            %Grab the trim data
            ID0      = obj.BulkDataID0;
            TrimData = awi.methods.Nastran.getTrimData(FEM, Aircraft, LoadCases, MassCases, ID0, obj.RefNode);
            
            GustData = obj.getGustData(TrimData,FlightPoints,LoadCases);
            
            %% Set up the directory and print model bulk data
            
            %Create a folder for the output
            if nargin < 4
                dirName = [datestr(now, 30), '_gust_analysis'];
                [bSuccess, ~, ~] = mkdir(pwd, dirName);
                assert(bSuccess == true, ['Unable to create output ', ...
                    'directory for the trim analsysis.']);
                outDir = fullfile(pwd, dirName);
            end
            
            %Export the FEM
            [~, includeFiles] = export(FEM, outDir, 'WriteHeaderFile', false);
            
            %% Write the header file
            gustFile = p.Results.DatFilename;
            if ~isempty(gustFile)
                gustFile =  obj.parseFilename(gustFile);
            elseif numel(LoadCases) > 1 || isempty(LoadCases.Name)
                gustFile = 'gust_analysis';
            else
                gustFile = obj.parseFilename(LoadCases.Name);
            end
            gustFile = fullfile(outDir, [gustFile, '.dat']);
            
            fid = fopen(gustFile, 'w');
            
            %Executive Control
            awi.fe.FEBaseClass.writeHeading(fid, 'E X E C U T I V E  C O N T R O L');
            fprintf(fid, 'SOL %i\r\n', 146);
            fprintf(fid, 'ECHOOFF         $ SUPPRESSES THE ECHO OF EXECUTIVE CONTROL\r\n');
            fprintf(fid, 'CEND\r\n');

            %% SPC temp ---------------------------------------
            fprintf(fid, 'SPC  = %i\r\n', ID0);
            % -------------------------------------------------

            %Case Control
            awi.fe.FEBaseClass.writeHeading(fid, 'C A S E  C O N T R O L');
            awi.fe.FEBaseClass.writeSubHeading(fid, 'O U T P U T  O P T I O N S');
            fprintf(fid, 'LINE = 99999   $ SPECIFIES THE NUMBER OF LINES PER PRINTED PAGE\r\n');
            fprintf(fid, 'ECHO = NONE    $ SUPPRESSES THE ECHO OF BULK DATA\r\n');
            %
            awi.fe.FEBaseClass.writeSubHeading(fid, 'O U T P U T  Q U A N T I T I E S');
            fprintf(fid, 'DISP(PLOT)  = ALL $ OUTPUT ALL DISPLACEMENTS\r\n');
            fprintf(fid, 'DISP = ALL $ OUTPUT ALL DISPLACEMENTS\r\n');
            fprintf(fid, 'FORCE(PLOT) = ALL $ OUTPUT ALL ELEMENT FORCES\r\n');
            fprintf(fid, 'AEROF       = ALL $ OUTPUT ALL AERO FORCES\r\n');
            fprintf(fid, 'APRESSURE   = ALL $ OUTPUT ALL AERO PRESSURES\r\n');
            %
            
            fprintf(fid, 'METHOD = %i\r\n', GustData.EIGRL.SID );
            fprintf(fid, 'SDAMP  = %i\r\n', GustData.TABDMP1.TID);
            fprintf(fid, 'FREQ   = %i\r\n', GustData.FREQ.id);
            
            num=height(GustData.bulk_data);
            
            for caseID=1:num
                awi.fe.FEBaseClass.writeSubHeading(fid, 'S U B C A S E S')
                fprintf(fid, 'SUBCASE %i\r\n', caseID);
                fprintf(fid, 'LABEL = %s\r\n','Gust response');
                fprintf(fid, 'TSTEP  = %i\r\n', GustData.TSTEP.id);
                fprintf(fid, 'GUST  = %i\r\n', GustData.GUSTid+(caseID-1)*10);
                
            end
            
            awi.fe.FEBaseClass.writeSubHeading(fid, 'G L O B A L  C A R D S')
            
            %Bulk Data
            obj.defaultBulkStatement(fid);
            fprintf(fid, 'PARAM,MACH,%8.2f\r\n',FlightPoints.Mach);
            fprintf(fid, 'PARAM,Q,%8.2f\r\n',FlightPoints.DynPressure);
            fprintf(fid, 'PARAM,GUSTAERO,-1\r\n');
            
            %Write the cards
            awi.fe.FEBaseClass.writeColumnDelimiter(fid, 'normal');

             % SPC 
             fprintf(fid, '%-8s%-8i%-8i%-8i\r\n', 'SPC1', ID0, 123456, 1068);
            
             % frequency 
             fprintf(fid, '%-8s%-8i%-8.3f%-8.3f%-8i\r\n', 'FREQ1',GustData.FREQ.id,GustData.FREQ.F1,GustData.FREQ.DF,GustData.FREQ.NDF);
             
             % TSTEP
             fprintf(fid, '%-8s%-8i%-8i%-8.2f%-8i\r\n', 'TSTEP',GustData.TSTEP.id,GustData.TSTEP.N,GustData.TSTEP.DT,GustData.TSTEP.NO1);
             
             % TABDMP
             fprintf(fid, ['%-8s%-8i%-8s\r\n%-8s', repmat('%-8s', [1, numel(GustData.str)]), '\r\n'], ...
                 'TABDMP1', GustData.TABDMP1.TID, 'CRIT', blanks(8), GustData.str{:});
                       
             % DAREA
             fprintf(fid, '%-8s%-8i%-8i%-8i%-8.2f\r\n', 'DAREA',GustData.DAREA.sid,GustData.DAREA.p1,GustData.DAREA.c1,GustData.DAREA.a1);
             
             awi.fe.FEBaseClass.writeColumnDelimiter(fid, 'normal');
             
              % Gust definitions for each subcase 
              for caseID=1:num
                  
                  % TLOAD1
                  fprintf(fid, ['%-8s%-8i%-8i',blanks(16),'%-8i\r\n'], 'TLOAD1',GustData.TLid+caseID-1,GustData.DAREA.sid,GustData.TABLED1.id+caseID-1);
                  
                  % - DLOAD
                  %   fprintf(fid, '%-8s%-8i%-8.3f%-8.3f%-8i\r\n', 'DLOAD',gust_data.DLid,gust_data.DLS,gust_data.DLSi,gust_data.TLid);
                  
                  % Gust card
                  fprintf(fid, '%-8s%-8i%-8i%-8.3f%-8.1f%-8.1f\r\n', 'GUST',GustData.GUSTid+(caseID-1)*10,...
                      GustData.TLid+caseID-1,GustData.Wg(caseID),GustData.X0,FlightPoints.AcVelocity);
                  
                  % TABLED1
                  fprintf(fid, '%-8s%-8i\r\n', 'TABLED1',GustData.TABLED1.id+caseID-1);
                  
                  fprintf(fid, '        %-8.3f%-8.3f%-8.3f%-8.3f%-8.3f%-8.3f%-8.3f%-8.3f\r\n',GustData.bulk_data(caseID,:));
                  
                  fprintf(fid, [blanks(8),'%-8s\r\n'],'ENDT');
                  
                  awi.fe.FEBaseClass.writeColumnDelimiter(fid, 'normal');
                  
              end
              
      
             %   EIGR             
             fprintf(fid, '%-8s%-8i%-8s%-#8.3g%-16s%-8i\r\n', 'EIGR', GustData.EIGRL.SID, 'MGIV', 0, blanks(16),GustData.EIGRL.ND);
             
%              fprintf(fid, '%-8s%-8i%-8s%-#8.3g%-#8.3g\r\n', 'EIGR', GustData.EIGRL.SID, 'MGIV', 0, 1);
             
             % - AERO
             obj.writeUnsteadyAeroEntry(fid, Aircraft, FlightPoints(1));
             
             
             % - MKAERO1
             GustData.MKAERO.M = unique([FlightPoints.Mach]);
             GustData.MKAERO.K = [0.001, 0.005, 0.01, 0.03, 0.06, 0.1, 0.2, 0.4, 0.7, 1.0, 1.3, 1.6, 2, 2.5, 3, 3.5];
             obj.writeMKAERO1(fid, GustData.MKAERO);
                        
            %Additional bulk data files
            if ~iscell(includeFiles)
                includeFiles = {includeFiles};
            end
            awi.fe.FEBaseClass.writeSubHeading(fid, 'I N C L U D E  F I L E S');
            awi.fe.FEBaseClass.writeIncludeStatement(fid, includeFiles);
            
            %End of file
            fprintf(fid, 'ENDDATA\r\n');
            
            %Close the file
            fclose(fid);                        
            
        end
    end
    
    methods % flutter analysis
        function FlutterResults = flutter(obj, Aircraft, FlightPoints, MassCases, varargin)
            %flutter Executes a Solution 145 analysis in MSC.Nastran using the
            %analysis model 'obj.AnalysisModel', the aircraft data in
            %'Aircraft', the flight point data in 'FlightPoints' and any
            %additional mass data in 'MassCases'.
            
            %Parse
            validateattributes(Aircraft    , {'awi.model.Aircraft'}   , ...
                {'scalar'}, 'flutter', 'Aircraft');
            validateattributes(FlightPoints, {'awi.model.FlightPoint'}, ...
                {'row'}   , 'flutter', 'FlightPoints');
            if nargin < 4
                MassCases = [];
            else
                validateattributes(MassCases, {'awi.model.MassCase'}, ...
                    {'row'}, 'flutter', 'MassCases');
            end
            p = inputParser;
            addParameter(p, 'BoundaryCondition', '123456', @(x)any(validatestring(x, {'123456', '1246'})));
            parse(p, varargin{:});
            
            FlutterResults = [];
            
            %Write the header file
            flutterFile = obj.writeFlutterFile(Aircraft, FlightPoints, MassCases, p.Results.BoundaryCondition);
            
            %Run the analysis
            %   - Need to change the current directory to force Nastran to
            %     output the files in the correct location
            %             obj.runNastran(flutterFile);
            
            %Extract the results from the Nastran output
            %   - The aero forces & pressures are in the .f06
            %             NastranResults = obj.extractNastranResults(flutterFile, 'ReadF06', true);
            
            %Generate the 'awi.model.TrimResult' object
            %             TrimResults = obj.generateAWIResults(NastranResults);
            
        end
        function flutterFile = writeFlutterFile(obj, Aircraft, FlightPoints, MassCases, spc, outDir, varargin)
            %writeFlutterFile Writes a header file for a aeroelastic flutter
            %analysis using SOL 145 for the flight point data defined
            %'FlightPoints' using any additional mass data in 'MassCases'.
                        
            %Parse
            if nargin < 4
                MassCases = [];
            end
            if nargin < 5 || isempty(spc)
                spc = '123456';
            else
                validatestring(spc, {'123456', '1246'});
            end            
            p = inputParser; 
            addParameter(p, 'RequestModeshapes', false , @(x)validateattributes(x, {'logical'}, {'scalar'}));
            addParameter(p, 'FlutterMethod'    , 'PKNL', @(x)any(validatestring(lower(x), {'pk', 'pknl'})));
            addParameter(p, 'DatFilename'      , []    , @(x)validateattributes(x, {'char'}, {'row'}));
            parse(p, varargin{:});
            
            %Create a folder for the output
            if nargin < 6 || isempty(outDir)
                dirName = [datestr(now, 30), '_flutter_analysis'];
                [bSuccess, ~, ~] = mkdir(pwd, dirName);
                assert(bSuccess == true, ['Unable to create output ', ...
                    'directory for the flutter analsysis.']);
                outDir = fullfile(pwd, dirName);
            else
                assert(isfolder(outDir), ['Expected the output ', ...
                    'directory to be a valid folder path'])
            end
            
            %Filename?
            flutterFile = p.Results.DatFilename;
            if isempty(flutterFile)
                flutterFile = 'flutter_analysis';
            else
                flutterFile =  obj.parseFilename(flutterFile);                
            end
                        
            %Grab the FEM
            FEM = obj.AnalysisModel;
            
            %Grab the flutter data
            RefNode_         = obj.RefNode;
            ID0              = obj.BulkDataID0;
            if isempty(RefNode_)
                %                 error(['Unable to proceeed as there is no reference ', ...
                %                     'node for the constraints to be applied to.']);
                RefNode_ = obj.WingNode(1);
            end
            BulkData.SPCdof  = spc;
            BulkData.RefGrid = RefNode_;
            FlutterData      = obj.getFlutterData(FlightPoints, BulkData, ID0);
            
            %Set velocities to negative if modeshapes are requested
            if p.Results.RequestModeshapes
                idx = ismember({FlutterData.FLFACT.LAB}, 'VELOCITY');
%                 FlutterData.FLFACT(idx).F = -FlutterData.FLFACT(idx).F;  
                velocity_input=unique([FlightPoints.AcVelocity,FlutterData.FLFACT(idx).F]);
                
                % only write modes of interest (at particular velocity)
                velocity_input(velocity_input==FlightPoints.AcVelocity)=-FlightPoints.AcVelocity;
                FlutterData.FLFACT(idx).F = velocity_input;
            end
            
            %Ensure unique (V,M,rho) combinations for PK method
            if strcmpi(p.Results.FlutterMethod, 'pk')
                for i = 1 : numel(FlutterData.FLFACT)
%                     FlutterData.FLFACT(i).F = unique(FlutterData.FLFACT(i).F);
                    FlutterData.FLFACT(i).F = FlutterData.FLFACT(i).F;
                end
            end
            
            %Use correct flutter method 
            FlutterData.FLUTTER.METHOD = upper(p.Results.FlutterMethod);
            
            %% Print model bulk data
            
            %Export the FEM
            [~, includeFiles] = export(FEM, outDir, 'WriteHeaderFile', false);
            
            %% Write the header file
            flutterFile = fullfile(outDir, [flutterFile, '.dat']);
            
            fid = fopen(flutterFile, 'w');
            
            %Executive Control
            awi.fe.FEBaseClass.writeHeading(fid, 'E X E C U T I V E  C O N T R O L');
            fprintf(fid, 'SOL %i\r\n', 145);
            fprintf(fid, 'ECHOOFF         $ SUPPRESSES THE ECHO OF EXECUTIVE CONTROL\r\n');
            fprintf(fid, 'CEND\r\n');
            
            %Case Control
            awi.fe.FEBaseClass.writeHeading(fid, 'C A S E  C O N T R O L');
            awi.fe.FEBaseClass.writeSubHeading(fid, 'O U T P U T  O P T I O N S');
            fprintf(fid, 'LINE = 99999999   $ SPECIFIES THE NUMBER OF LINES PER PRINTED PAGE\r\n');
            fprintf(fid, 'ECHO = NONE       $ SUPPRESSES THE ECHO OF BULK DATA\r\n');
            %
            awi.fe.FEBaseClass.writeSubHeading(fid, 'O U T P U T  Q U A N T I T I E S');
            if p.Results.RequestModeshapes
                fprintf(fid, 'DISP(PLOT)  = ALL $ OUTPUT ALL DISPLACEMENTS\r\n');
                fprintf(fid, 'DISP = ALL $ OUTPUT ALL DISPLACEMENTS\r\n');
                fprintf(fid, 'SVEC = ALL $ OUTPUT ALL DISPLACEMENTS\r\n');
                fprintf(fid, 'VECTOR(SORT1,REAL)  = ALL $ OUTPUT ALL EIGEN VECTORS\r\n');
            end
            %             fprintf(fid, 'FORCE(PLOT) = ALL $ OUTPUT ALL ELEMENT FORCES\r\n');
            %             fprintf(fid, 'AEROF       = ALL $ OUTPUT ALL AERO FORCES\r\n');
            %             fprintf(fid, 'APRESSURE   = ALL $ OUTPUT ALL AERO PRESSURES\r\n');
            %
            awi.fe.FEBaseClass.writeSubHeading(fid, 'G L O B A L  C A R D S')
            fprintf(fid, 'SPC  = %i\r\n', FlutterData.SPC.SID);
            %
            awi.fe.FEBaseClass.writeSubHeading(fid, 'S U B C A S E S')
            pad   = blanks(2);
            scNo  = 1;
            bDisp = false;
            obj.printFlutterSubcases(fid, FlutterData, [], scNo, pad, bDisp);
            
            %Bulk Data
            obj.defaultBulkStatement(fid);
            %   - Flutter data
            obj.writeFlutterEntries(fid, FlutterData);
            %   - AERO
            obj.writeUnsteadyAeroEntry(fid, Aircraft, FlightPoints(1));
%             %charles added
%             obj.writeSteadyAeroEntry(fid, Aircraft);
            
            % remove any SUPORT or SPC for free - free condition (charles)
%             %   - SUPORT
%             if strcmp(FlutterData.SPC.C, '123456')
%                 fprintf(fid, '%-8s%-8i%-8i%-8i\r\n', 'SPC1', ID0, 123456, FlutterData.SPC.G);
% %                 fprintf(fid, '%-8s%-8i%-8s\r\n', 'SUPORT', FlutterData.SPC.G(34).GID, '35');
%             end

            %   - SPC
            mni.printing.cards.SPC1(ID0, str2num(spc), [FlutterData.SPC.G]).writeToFile(fid);


            %   - Additional bulk data files
            if ~iscell(includeFiles)
                includeFiles = {includeFiles};
            end
            awi.fe.FEBaseClass.writeSubHeading(fid, 'I N C L U D E  F I L E S');
            awi.fe.FEBaseClass.writeIncludeStatement(fid, includeFiles);
            
            %End of file
            fprintf(fid, 'ENDDATA\r\n');
            
            %Close the file
            fclose(fid);
            
        end
    end
    
    methods % structural optimisation
        function OptResults = structuralOptimisation(obj, Aircraft, LoadCases, Options, optPartName, optPartCode, outdir)
            %structuralOptimisation Executes a Solution 200 analysis in
            %MSC.Nastran.
            %
            % The geometry model and loadcases are used to define the
            % Finite Element Model (FEM) and analysis subcases, typically
            % these are aeroelastic but flutter and global buckling can
            % also be considered.
            %
            % The 'Options' variable is an instance of
            % 'awi.methods.BeamOptParam' and defines the type of box-beam
            % idealisation that is used to represent the internal
            % structure of the lifting surfaces. It also controls various
            % parameters of the SOL 200 analysis.
            %
            % 'optPartName' and 'optPartCode' define the parts of the
            % geometry model that will be optimised.
            
            OptResults = [];
            
            bOpt = Nas.canOptimise;
            if ~bOpt %Escape route
                return
            end
            
            %Parse inputs
            
            %Write the optimisation file
            optFile = writeSol200File(obj, Aircraft, LoadCases, Options, optPartName, optPartCode, outdir);
            
            %Execute analysis
            %   - Need to change the current directory to force Nastran to
            %     output the files in the correct location
            obj.runNastran(optFile);
            
            %Generate results object
            
        end
        function [optFile, varargout] = writeSol200File(obj, Aircraft, LoadCases, Options, optPartName, optPartCode, outdir, varargin)
            %writewriteSol200File
            
            %Parse inputs
            stck = dbstack;
            funcName = strsplit(stck(1).name, '.');
            funcName = funcName{end};
            if nargin < 6
                inputs = {'Aircraft', 'LoadCases', 'Options', 'optPart', 'optCode'};
                nIn    = numel(inputs);
                msg    = ['The function ''%s'' requires %i inputs:\n\n', ...
                    repmat('\t%s\n', [1, nIn])];
                error(msg, funcName, nIn, inputs{:});
            end
            if nargin < 7 || isempty(outdir)
                outdir = pwd;
            end
            p = inputParser;
            addRequired(p , 'Aircraft'   , @(x)validateattributes(x, {'awi.model.Aircraft'}      , {'scalar'}, funcName, 'Aircraft'));
            addRequired(p , 'LoadCases'  , @(x)validateattributes(x, {'awi.model.LoadCase'}      , {'vector'}, funcName, 'LoadCases'));
            addRequired(p , 'Options'    , @(x)validateattributes(x, {'awi.methods.BeamOptParam'}, {'scalar'}, funcName, 'Options'));
            addRequired(p , 'optPartName', @iscellstr);
            addRequired(p , 'optPartCode', @iscellstr);
            addRequired(p , 'outdir'     , @isfolder);
            addParameter(p, 'WriteHeaderFile', true, @(x)validateattributes(x, {'logical'}, {'scalar'}, funcName, 'WriteHeaderFile'));
            addParameter(p, 'WriteDataFiles' , true, @(x)validateattributes(x, {'logical'}, {'scalar'}, funcName, 'WriteDataFiles'));
            parse(p, Aircraft, LoadCases, Options, optPartName, optPartCode, outdir, varargin{:});
            
            %% Determine whether the Aircraft can be optimised
            %   - Material properties defined?
            %   - Cross-sections defined?
            
            if Options.UseBeamCrossSection
                error(['SOL 200 does not return correct stress when ', ...
                    'using PBEAML. Use PBEAM instead. (User error - ', ...
                    'Not NASTRAN error).']);
            end
            
            %% Prepare geometry model for conversion to FEM
            %   - Apply symmetry conditions
            %   - Update rib/ stringer pitch
            %   - Apply global parameters (UseBeamCrossSection)
            %   - Generate 'awi.model.BoxBeam' objects
            %   - Update dimensions of box-beam objects based on the
            %     prescribed geometrical cross-section of each
            %     'awi.model.Beam' object.
            %   - Convert to FE model
            %   - Remove material density where required
            
            if Options.ModelXZSymmetry
                [optPartName, optPartCode] = i_modelSymmetry(Aircraft, optPartName, optPartCode);
            end
            
            function [optPartName, optPartCode] = i_modelSymmetry(Aircraft, optPartName, optPartCode)
                %i_modelSymmetry Removes any port-side objects from the
                %model and sets 'ModelAero' to false for the VTP element.
                %(i.e. Removes any aerodynamic panels on the VTP from the
                %FEM)
                %
                % TODO - Relies on specific naming conventions for the
                % parts. Improve this!
                
                parts     = flatlist(Aircraft);
                pNames = {parts.Name};
                
                idxPort   = contains(pNames, 'Port');
                delete(parts(idxPort));
                parts = parts(~idxPort);
                
                idxPort = contains(optPartName, 'Port');
                optPartName(idxPort) = [];
                optPartCode(idxPort) = [];
                
                vtp = parts(arrayfun(@(o) and(isa(o, 'awi.model.LiftingSurface'), strcmp(o.Name, 'VTP')), parts));
                vtp.ModelAero = false;
                
            end
            
            %TODO - Remove this and make it an input to the method
            OptPart = i_getOptPart(Aircraft, optPartName);
            if isempty(OptPart)
                warning('Unable to find any of the ''awi.model'' objects using the part names provided');
                return
            end
            
            function OptPart = i_getOptPart(Aircraft, optPartName)
                %i_getOptPart Returns the handle to the 'awi.model' objects
                %that will be included in the Design Model.
                %
                % TODO - Currently only valid for lifting surfaces.
                
                parts   = flatlist(Aircraft);
                LS      = parts(arrayfun(@(o) isa(o, 'awi.model.LiftingSurface'), parts));
                OptPart = LS(ismember({LS.Name}, optPartName));
                
            end
            
            %Update box-beam parameters based on input options
            prp = {'RibPitch', 'StringerPitch', 'UseBeamCrossSection'};
            val = get(Options, prp);
            set(OptPart, prp, val);
            
            %Generate BoxBeam objects along the cross-section
            [thickVar, areaVar, optVarNames] = i_grabBoxComponentNames(Options);
            arrayfun(@(o) createBoxBeamObjects(o, Options.BoxBeamType, ...
                thickVar, Options.SparHeight), OptPart);
            
            function [thickVar, areaVar, optVarNames] = i_grabBoxComponentNames(Options)
                %i_grabBoxComponentNames Returns the names of the box
                %parameters that will be optimised as well as the value of
                %the thickness and area dimensions.
                
                areaVar = {};
                
                bb = awi.model.BoxBeam;
                
                %Grab box properties relating to the current box type
                optVarNames = bb.CurrentPropNames;
                
                %The width & height are fixed during the optimisation
                optVarNames(ismember(optVarNames, {'Width', 'Height'})) = [];
                
                %Initial thickness & area for a single section
                thickVar = repmat(Options.InitialThick, [1, numel(optVarNames)]);
                
            end
            
            %Convert to FEM
            FEM       = convertToFE(Aircraft);
            allFEM    = flatlist(FEM);
            
            %Create unique mass model
            [MassCases, ~, ind] = unique([LoadCases.MassCase], 'stable');
            if isempty(MassCases)
                MassModel = [];
                MassSets  = [];
            else
                %Convert to FE
                MassModel = convertThisToFE(MassCases, allFEM);
                %Create a 'awi.fe.MassSet' object for each loadcase
                MassSets = arrayfun(@(~) awi.fe.MassSet, 1 : numel(ind));
                for iM = 1 : numel(MassCases)
                    set(MassSets(ind == iM), 'MassModel', MassModel(iM));
                end
                %Update ID numbers now that extra masses have been generated
                assignIDnumbers([allFEM, MassModel], MassSets);
            end
            
            %Retain the FEM relating to the optimisation
            OptFEM = allFEM(ismember([allFEM.GeometryObject], OptPart));
            
            %Strip material density from all 'awi.fe.Material cards'
            %related to the optimisation model.
            %   - Element mass is given by beam NSM instead.
            Mat    = [OptFEM.Materials];
            set(Mat, 'Rho', 0);
            
            obj.AnalysisModel = FEM;
            
            %% Define analysis data
            
            [TrimData, StaticData, BuckData, FlutterData] = obj.getAnalysisData(FEM, Aircraft, LoadCases, MassSets, Options);
            
            %% Define the objective function
            
            %Total model mass (excluding residual mass models)
            ModelMass = awi.fe.opt.DesignResponse;
            ModelMass.LABEL = 'WGHT';
            ModelMass.ResponseType = 'WEIGHT';
            ModelMass.ResponseAttributes = {3, 3, []};
            
            WeightResp = ModelMass;
            
            %             %Define equation normalising the model mass by the SUGAR AUM
            %             %   - TODO : Remove this hardcoded bit and replace with a
            %             %   variable.
            %             M = 39598.58;
            %             MassEqn     = awi.fe.opt.DesignEquation;
            %             MassEqn.Eqn = ['JX(MASS) = MASS / ', num2str(M, '%.4f')];
            %
            %             %Objective function
            %             WeightResp       = awi.fe.opt.DesignResponseEqn;
            %             WeightResp.LABEL = 'JX';
            %             WeightResp.DesignEquation  = MassEqn;
            %             WeightResp.DesignResponses = ModelMass;
            %             addFEData(FEM(1), WeightResp); %, MassEqn, ModelMass);
            
            addFEData(FEM(1), WeightResp); %, MassEqn, ModelMass);
            
            %% Define the design equations
            %   - Box geoemtrical properties (Ixx, Izz, A, J, etc.)
            %   - Box internal area
            %   - Direct and shear stresses
            %   - Von Mises stress
            %   - Reserve factors
            
            %Design Equations
            %   - Define equations for the beam properties based on the
            %     specified cross-section
            [BoxEqn, BoxArea, SkinThick, StrArea, StrThick, StrZna, StrMomArea] = i_defineBoxPropEqn(Options.BoxBeamType);
            %   - Define equations for beam stresses based on specified cross-section
            BoxStressEqn = i_defineBoxStressEqn(Options.BoxBeamType);
            %   - Define an equation that takes in the direct and shear
            %     stress and returns the Von Mises Stress
            VMEqn     = awi.fe.opt.DesignEquation;
            VMEqn.Eqn = 'VM(SIG, TAU) = SQRT(SIG**2 + 3.0*TAU**2)';
            %   - Define an equation that takes in a stress value and
            %     divides it by the ultimate stress to calculate the
            %     Resrve Factor
            RFEqn     = awi.fe.opt.DesignEquation;
            RFEqn.Eqn = 'RF(ULT, ACTUAL) = (ULT/(ACTUAL + 0.001))';
            %   - Define an equation that expresses the panel instability
            PanBuckEqn = awi.fe.opt.DesignEquation;
            PanBuckEqn.Eqn = 'SIGBUCK(PI, KC, BSTR, E, NU, TSK) = ((PI**2 * E * KC) / (12.0 * (1.0 - NU**2))) * (TSK / BSTR)**2';
            %   - Define an equation that expresses the stringer instability
            StrBuckEqn = awi.fe.opt.DesignEquation;
            StrBuckEqn.Eqn = 'SIGBUCK(PI, E, LR, ASTR, ISTR) = (PI**2 * E * ISTR) / (LR**2 * ASTR)';
            
            %Always need the box equations
            addFEData(FEM(1), BoxEqn);
            
            %Only require the other equations if SOL200 is doing all the
            %optimisation
            if ~Options.ExitAfterSensitivity
                addFEData(FEM(1), BoxArea, SkinThick, StrArea, ...
                    StrThick, StrZna, StrMomArea, BoxStressEqn, VMEqn, RFEqn, ...
                    PanBuckEqn, StrBuckEqn);
            end
            
            function [BoxEqn, BoxArea, SkinThick, StrArea, StrThick, StrZna, StrMomArea] = i_defineBoxPropEqn(boxType)
                %i_defineBoxPropEqn Defines equations relating the box
                %dimensions to the box properties.
                %
                % The number of equations provided will be determined by
                % the type of box.
                
                switch boxType
                    case 'SymmetricBox'
                        %Cross-Sectional Area (A)
                        AreaEqn     = awi.fe.opt.DesignEquation;
                        AreaEqn.Eqn = 'Wi(TC, TS, W, H) = W - (2.0 * TS);Hi = H - (2.0 * TC);A = (W * H) - (Wi * Hi)';
                        %Second Moment of Area (Izz)
                        IzzEqn     = awi.fe.opt.DesignEquation;
                        IzzEqn.Eqn = 'Wi(TC, TS, W, H) = W - (2.0 * TS);Hi = H - (2.0 * TC);Izz = ((H * W**3) / 12.0) - ((Hi * Wi**3) / 12.0)';
                        %Second Moment of Area (Ixx)
                        IxxEqn     = awi.fe.opt.DesignEquation;
                        IxxEqn.Eqn = 'Wi(TC, TS, W, H) = W - (2.0 * TS);Hi = H - (2.0 * TC);Ixx = ((W * H**3) / 12.0) - ((Wi * Hi**3) / 12.0)';
                        %Polar Moment of Area (J)
                        JEqn     = awi.fe.opt.DesignEquation;
                        JEqn.Eqn = 'J(TC, TS, W, H) = ((2.0 * TC * TS) * (W - TS)**2 * (H - TC)**2) / (W * TS + H * TC - TC**2 - TS**2)';
                        %Non-structural Mass (NSM)
                        NSMEqn     = awi.fe.opt.DesignEquation;
                        NSMEqn.Eqn = 'Wi(TC, TS, W, H, RHO) = W - (2.0 * TS);Hi = H - (2.0 * TC);NSM = RHO * (W * H - Wi * Hi)';
                        %Non-structural Inertia (NSI)
                        NSIEqn     = awi.fe.opt.DesignEquation;
                        NSIEqn.Eqn = ['Wi(TC, TS, W, H, RHO) = W - (2.0 * TS);Hi = H - (2.0 * TC);MC = RHO * TC * Wi;RC = (H-TC)/2.0;', ...
                            'IC = (MC/12.0)*(Wi**2 + tC**2) + MC*RC**2;MS = RHO * TS * H;RS = (W-TS)/2.0;IS = (MS/12.0)*(H**2 + TS**2) + MS*RS**2;', ...
                            'I = 2.0 * (IC + IS)'];
                        %Collect
                        BoxEqn = [AreaEqn, IzzEqn, IxxEqn, JEqn, NSMEqn, NSIEqn];   %Order of eqns must match order of properties in 'beamProps'
                        %Area enclosed by the box centre-line
                        BoxArea = awi.fe.opt.DesignEquation;
                        BoxArea.Eqn = 'Ai(TC, TS, W, H) = (W - TS) * (H - TC)';
                        
                        %case 'Box'
                        %xNAEqn =
                        %zNAEqn =
                    otherwise
                        error('Unknown Box Beam Type provided - Unable to run optimisation.');
                end
                
                %Thickness of the skin
                SkinThick     = awi.fe.opt.DesignEquation;
                SkinThick.Eqn = 'TSKN(TC, ASBT) = TC / (ASBT + 1.0)';
                
                %Area of stringer
                StrArea     = awi.fe.opt.DesignEquation;
                StrArea.Eqn = 'ASTR(TC, BSTR, ASBT) = (ASBT * BSTR * TC) / (ASBT + 1.0)';
                
                %Thickness of stringer web
                StrThick = awi.fe.opt.DesignEquation;
                StrThick.Eqn = 'TSTR(WSTR, ASTR, TSKN) = ((2.0 * ASTR) - (WSTR * TSKN)) / (6.0 * TSKN)';
                
                %Stringer neutral axis position
                StrZna = awi.fe.opt.DesignEquation;
                StrZna.Eqn = ['AWEB(WSTR, TSKN, TSTR) = TSTR * 3.0 * TSKN; ' , ...
                    'ZWEB = 2.0 * TSKN; AFLANGE = 0.5 * TSKN * WSTR; ', ...
                    'ZFLANGE = TSKN / 4.0; STRZNA = ((AWEB * ZWEB) + '        , ...
                    '(AFLANGE * ZFLANGE)) / (AWEB + AFLANGE)'];
                
                %Second moment of area of the stringer (Izz)
                StrMomArea = awi.fe.opt.DesignEquation;
                StrMomArea.Eqn = ['AWEB(WSTR, TSKN, TSTR, STRZNA) = '    , ...
                    'TSTR * 3.0 * TSKN; ZWEB = 2.0 * TSKN; '     , ...
                    'AFLANGE = 0.5 * TSKN * WSTR; ZFLANGE = TSKN / 4.0; ', ...
                    'ISTR = ((TSKN * (3.0 * TSTR)**3) / 12.0) + (AWEB * ', ...
                    '(ZWEB - STRZNA)**2) + (AFLANGE * '                  , ...
                    '(ZFLANGE - STRZNA)**2) + (WSTR * (0.5 * TSKN)**3 / 12.0)'];
                
            end
            
            function BoxStressEqn = i_defineBoxStressEqn(boxType)
                %i_defineBoxPropEqn Defines equations relating the box
                %dimensions to the box stress distribution.
                %
                %The number of equations provided will be determined by the
                %type of box.
                
                switch boxType
                    case 'SymmetricBox'
                        %Direct stress
                        DStressEqn     = awi.fe.opt.DesignEquation;
                        DStressEqn.Eqn = 'SIG(TS, TC, W, H, M1, M2, AF, A, I1, I2) = ((ABS(M1)/I1)*((W/2.0) - (TS/2.0))) + ((ABS(M2)/I2)*((H/2.0) - (TC/2.0))) + (AF/A)';
                        %Shear Stress
                        SStressEqn     = awi.fe.opt.DesignEquation;
                        SStressEqn.Eqn = 'TAU(T, W, H, K, F1, F2, TQ, Ai) = (ABS(TQ)/(2.0*Ai*T)) + (ABS(F1)/(2.0*H*T)) + (K*(ABS(F2)/(2.0*W*T)))';
                        %Collect
                        BoxStressEqn = [DStressEqn, SStressEqn];
                end
                
            end
            
            %% Define constants (box widths, material strengths, etc.)
            
            i_defineMatProps(OptFEM, optPartCode);
            
            function MatTable = i_defineMatProps(OptFEM, optPartCode)
                %i_defineMatProps Defines the 'awi.fe.opt.DesignTable'
                %objects containing information aboutthe material
                %properties for each part that is being optimised using the
                %information stored in the 'awi.model.Material' object
                %attached to the FEM geometry object.
                %
                % Properties defined:
                %   - Ultimate tensile strength (STE)
                %   - Ultimate shear strength (SSH)
                %   - Material density (RHO)
                
                for i = 1 : numel(OptFEM)
                    MatTable = awi.fe.opt.DesignTable;
                    AWIMat   = OptFEM(i).GeometryObject.Material(1);
                    MatTable.Prop = strcat({'STE_', 'SSH_', 'RHO_', 'E_', 'NU_'}, optPartCode{i});
                    MatTable.Val  = cellstr(num2str( ...
                        [AWIMat.Sigma11_T, AWIMat.Tau12, AWIMat.Rho, AWIMat.E, AWIMat.Nu]', '%-#.3g'))';
                    %^^ Specify the string directly to ensure we fit into 8
                    %character limit
                    addFEData(OptFEM(i), MatTable);
                end
                
            end
            
            BoxDims = i_defineBoxHeightAndWidth(OptFEM, optPartCode);
            
            function BoxDims = i_defineBoxHeightAndWidth(OptFEM, optPartCode)
                %i_defineBoxHeightAndWidth Defines the
                %'awi.fe.opt.DesignTable' objects containing information
                %about the height and width of each box section.
                
                BoxDims = arrayfun(@(~) awi.fe.opt.DesignTable, 1 : numel(OptFEM));
                
                for i = 1 : numel(OptFEM)
                    
                    %Calculate eta position of beams
                    beamEta = calculateBeamEta(OptFEM(i));
                    beamEta = [beamEta(1, :), beamEta(2, end)];
                    
                    %Interpolate the height & width to beam positions
                    bb_eta  = OptFEM(i).GeometryObject.BoxBeam_eta;
                    height  = [OptFEM(i).GeometryObject.BoxBeam.Height];
                    width   = [OptFEM(i).GeometryObject.BoxBeam.Width];
                    h       = interp1(bb_eta, height, beamEta);
                    w       = interp1(bb_eta, width , beamEta);
                    
                    %Construct property codes
                    hCode   = arrayfun(@(n) sprintf('%s_H%i', optPartCode{i}, n), 1 : numel(h), 'Unif', false);
                    wCode   = arrayfun(@(n) sprintf('%s_W%i', optPartCode{i}, n), 1 : numel(w), 'Unif', false);
                    
                    %Assign to object
                    BoxDims(i).Prop = [hCode, wCode];
                    BoxDims(i).Val  = [h, w];
                    
                    addFEData(OptFEM(i), BoxDims(i));
                    
                end
                
            end
            
            %Additional constants
            %   - "KLOAD"  : Transverse load factor
            %   - "AS_BT"  : Stringer-to-skin area ratio
            %   - "B_STR"  : Stringer pitch
            %   - "L_RIB"  : Rib pitch
            %   - "PI"     : 3.14...
            %   - "KC"     : Buckling factor for panels
            %   - "KWING"  : Mass factor for the wing
            %   - "KSTRUT" : Mass factor for the strut
            %   - "W_STR"  : Width of the stringer base
            %   - "T_STR"  : Thickness of the stringer web
            OptConst      = awi.fe.opt.DesignTable;
            OptConst.Prop = {'KLOAD', 'AS_BT', 'B_STR', 'L_RIB', 'PI', 'KC', 'KWING', 'KSTRUT', 'W_STR'};
            OptConst.Val  = { ...
                num2str(Options.TransverseLoadFactor), ...
                num2str(Options.As_bT)               , ...
                num2str(Options.StringerPitch)       , ...
                num2str(Options.RibPitch)            , ...
                '3.14159', '4.0', '1.05', '1.3899', '0.0889'};
            
            addFEData(FEM(1), OptConst);
            
            %% Define the design variables
            
            %How many design variables in each box-beam section?
            nOptVar = numel(optVarNames);
            
            %How many beam properties in each FEM?
            nBP  = arrayfun(@(fm) numel(fm.BeamProps), OptFEM);
            
            %Using tapered beams with same cross-section at end-B_{i} and
            %end-A{i+1} so need (nBP+1) sets of design variables.
            nDV = (nBP + 1) .* nOptVar;
            ub  = cumsum(nDV);
            lb  = [1, ub(1 : end - 1)  + 1];
            DesVar = arrayfun(@(~) awi.fe.opt.DesignVariable, 1 : sum(nDV));
            
            %How are the design variables labelled?
            prpName  = OptPart(1).BoxBeam(1).CurrentPropNames;
            optLabel = OptPart(1).BoxBeam(1).CurrentPropLabel;
            optLabel(ismember(prpName, {'Width', 'Height'})) = [];
            
            for i = 1 : numel(OptFEM)
                
                %Assign design variables
                if Options.OptXZSymmetry && contains(optPartName{i}, 'Port')
                    stbdInd = returnStarboardIndex(optPartName, i);
                    dv      = DesVar(lb(stbdInd) : ub(stbdInd));
                else
                    dv  = DesVar(lb(i) : ub(i));
                    %Define the label for the design variables
                    lab = strcat(repmat(optPartCode(i), [nOptVar, 1]), optLabel');
                    lab = arrayfun(@(i) strcat(lab, {num2str(i)}), 1 : nBP(i) + 1, 'Unif', false);
                    lab = horzcat(lab{:});
                    set(dv, {'LABEL'}, lab(:));
                    %Assign bounds
                    dvlb = repmat(Options.MinThick, [numel(dv), 1]);
                    dvub = repmat(Options.MaxThick, [numel(dv), 1]);
                    %Grab box dimensions
                    w = str2double(BoxDims(i).Val(contains(BoxDims(i).Prop, [optPartCode{i}, '_W'])));
                    h = str2double(BoxDims(i).Val(contains(BoxDims(i).Prop, [optPartCode{i}, '_H'])));
                    dim(1 : 2 : numel(dv)) = h;
                    dim(2 : 2 : numel(dv)) = w;
                    %Check against the allowable box dimensions
                    idx = dvub > dim' ./ 2;
                    dvub(idx) = dim(idx) ./ 2;
                    %Set bounds
                    set(dv, {'XLB'}, num2cell(dvlb));
                    set(dv, {'XUB'}, num2cell(dvub));
                    addFEData(OptFEM(i), dv);
                    clear dim
                end
                
            end
            
            %Assign initial value
            %   - If any initial values exceed the upper bound then reset
            %   them to the upper bound.
            dv  = [OptFEM.DesVar];
            x0  = repmat(Options.InitialThick, [numel(dv), 1]);
            xub = vertcat(dv.XUB);
            idx = x0 >= xub;
            x0(idx) = xub(idx) .* 0.99; 
            set(dv, {'XINIT'}, num2cell(x0));
                                    
            %% Define the property relations, design responses & constraints for individual parts
            
            %Generate property codes for end-A & end-B
            if Options.UseBeamCrossSection
                beamProps = OptFEM(1).GeometryObject.BoxBeam(1).CurrentNasPropCode;
            else
                beamProps = OptFEM(1).GeometryObject.BoxBeam(1).CurrentNasPropNames;
            end
            nOptProp = numel(beamProps);
            pNam  = repmat(beamProps', [2, 1]);
            pNam  = strcat(pNam, [repmat({'(A)'}, [nOptProp, 1]) ; repmat({'(B)'}, [nOptProp, 1])]);
            nProp = numel(pNam);
            
            %Update shear centre locations
            i_updateShearCentreOffset(OptFEM);
            
            %Define property relations
            i_definePropRelation(OptFEM, optPartCode, optPartName, Options);
            
            function i_updateShearCentreOffset(OptFEM)
                %i_updateShearCentreOffset Updates the (x,z) offsets for
                %each FE beam element based on the location of the box-beam
                %shear centre.
                
                for i = 1 : numel(OptFEM)
                    %Grab 'awi.model.BoxBeam' objects and their positions
                    BB     = OptFEM(i).GeometryObject.BoxBeam;
                    bb_eta = OptFEM(i).GeometryObject.BoxBeam_eta;
                    eta    = calculateBeamEta(OptFEM(i));
                    eta    = [eta(1, :), eta(2, end)];
                    %Grab location of box-beam origin w.r.t the beam
                    org = [BB.SectionOrigin];
                    %Grab shear centre offsets (in local box-beam system)
                    xSC = [BB.xSC];
                    zSC = [BB.zSC];
                    %Shear centre offset in the global coordinate-system
                    off = org + [xSC ; zSC];
                    %As the port/starboard elements have the opposite
                    %reference for their beam coordinate system we have to
                    %flip the sign on the offset
                    if contains(OptFEM(i).Name, 'Starboard')
                        off(1, :) = abs(off(1, :));
                    end
                    %Interpolate to the FE beam positions
                    SCb = interp1(bb_eta, off', eta)';
                    %Account for end-A & end-B
                    SCx = [SCb(1, 1 : end - 1) ; SCb(1, 2 : end)];
                    SCz = [SCb(2, 1 : end - 1) ; SCb(2, 2 : end)];
                    %Assign to the 'awi.fe.Beam' objects
                    set(OptFEM(i).Beams, {'SCy'}, num2cell(SCx, 1)');
                    set(OptFEM(i).Beams, {'SCz'}, num2cell(SCz, 1)');
                end
                
            end
            
            function i_definePropRelation(OptFEM, optPartCode, optPartName, Options)
                %i_definePropRelation Defines the DVPREL2 entries that
                %define the box geometrical and inertial properties.
                
                for i = 1 : numel(OptFEM)
                    
                    nBP  = numel(OptFEM(i).BeamProps);
                    nBPr = nBP * nProp;
                    pr   = arrayfun(@(~) awi.fe.opt.PropRelationEqn  , 1 : nBPr);
                    
                    set(pr, 'TYPE', 'PBEAM');
                    
                    %Grab relevant 'awi.fe.opt' objects
                    if Options.OptXZSymmetry && contains(optPartName{i}, 'Port')
                        stbdInd = returnStarboardIndex(optPartName, i);
                        dv = OptFEM(stbdInd).DesVar;
                    else
                        dv   = OptFEM(i).DesVar;
                    end
                    [hLab, wLab] = i_grabBoxHeightAndWidth(OptFEM(i), optPartCode{i});
                    
                    %Assign handle to the beam property
                    bp = repmat(OptFEM(i).BeamProps, [nProp, 1]); %Duplicate handles for each end-A/end-B section
                    bp = bp(:);
                    set(pr, {'Property'}, num2cell(bp));
                    
                    %Assign the property name
                    pNam_ = repmat(pNam, [nBP, 1]);
                    set(pr, {'PropNameOrField'}, pNam_);
                    
                    %Assign the design variables
                    %   - Each section uses all design variables at that given section
                    dv_ = reshape(dv, [nOptVar, nBP + 1]);
                    %   - Assign each set of nOptVar design variables to the prop relation
                    dv_ = arrayfun(@(i) repmat(dv_(:, i), [1, nOptProp]), 1 : nBP + 1, 'Unif', false);
                    dv__ = cell(1, nBP * 2); %Account for End-A * End-B
                    dv__(1 : 2 : end - 1) = dv_(1 : end - 1);
                    dv__(2 : 2 : end)     = dv_(2 : end);
                    dv_ = horzcat(dv__{:});
                    set(pr, {'DesignVariables'}, num2cell(dv_, 1)');
                    
                    %Assign the DesOpt constants
                    dims  = [wLab ; hLab];
                    %   - Each section uses the width and height at the given section
                    const = cell(1, nBP * 2);
                    const(1 : 2 : end - 1) = arrayfun(@(i) dims(:, i)    , 1 : numel(wLab) - 1, 'Unif', false);
                    const(2 : 2 : end)     = arrayfun(@(i) dims(:, i + 1), 1 : numel(hLab) - 1, 'Unif', false);
                    %   - Account for all beam properties
                    const = arrayfun(@(i) repmat(const{i}, [1, nOptProp]), 1 : numel(const), 'Unif', false);
                    const = horzcat(const{:})';
                    set(pr, {'VarLabel'}, num2cell(const, 2));
                    %   - Update the 'VarLabel' for any NSM/NSI prop-relations so that the
                    %     material density is also passed to the eqn.
                    idx = or(contains(pNam_, 'NSM'), contains(pNam_, 'NSI'));
                    vl  = get(pr(idx), {'VarLabel'});
                    vl  = cellfun(@(x) [x, {['RHO_', optPartCode{i}]}], vl, 'Unif', false);
                    set(pr(idx), {'VarLabel'}, vl);
                    
                    %Assign the design equations
                    eqn = repmat(BoxEqn', [1, nBP * 2]);
                    set(pr, {'DesignEquation'}, num2cell(eqn(:)));
                    
                    addFEData(OptFEM(i), pr);
                    
                    
                    %Use this section if you want Nastran to output the
                    %value of the DVCREL cards as a DRESP2
                    %   - Allows you to check what Nastran is calculating
                    bMakeDRESP2 = false;
                    if bMakeDRESP2
                        %Make a DRESP2 entry for every DVPREL2 card
                        
                        DProp = pr(:);
                        
                        %Define new design responses based on these design properties
                        %   - Required because we cannot constrain a 'DVPREL' card
                        %   - Basically calculate the EI values again!
                        DRespProp = arrayfun(@(~) awi.fe.opt.DesignResponseEqn, 1 : numel(DProp));
                        lab   = strrep({DProp.PropNameOrField}, '(', '');
                        lab   = strrep(lab, ')', '');
                        index = repmat(1 : numel(lab) / nProp, [nProp, 1]);
                        index = index(:);
                        lab   = arrayfun(@(j) [optPartCode{i}, lab{j}, num2str(index(j), '%i')], 1 : numel(lab), 'Unif', false);
                        set(DRespProp, {'LABEL'}, lab');
                        set(DRespProp, {'DesignEquation'} , num2cell([DProp.DesignEquation]'));
                        set(DRespProp, {'DesignVariables'}, get(DProp, {'DesignVariables'}));
                        set(DRespProp, {'DesignConstants'}, get(DProp, {'VarLabel'}));
                        
                        %Make dummy constraints so these quantities are
                        %output
                        Constr = arrayfun(@(~) awi.fe.opt.DesignConstraint, 1 : numel(DRespProp));
                        set(Constr, {'DesignResponse'}, num2cell(DRespProp)');
                        set(Constr, 'LowerBound', -10);
                        
                        %Add it
                        addFEData(OptFEM(i), DRespProp, Constr)
                        
                    end
                    
                    
                end
                
            end
            
            %Define property responses & constraints for each part
            for iO = 1 : numel(OptFEM)
                
                %Check if we need to define responses for the port objects
                if Options.OptXZSymmetry && contains(optPartName{iO}, 'Port')
                    continue
                end
                
                nBP  = numel(OptFEM(iO).BeamProps);
                nBPr =  nBP * nProp;
                
                %Box dimensions
                [hLab, wLab] = i_grabBoxHeightAndWidth(OptFEM(iO), optPartCode{iO});
                dims  = [wLab ; hLab];
                
                %Material props
                %matProps = i_grabMaterialProps(OptFEM(iO),
                %optPartCode{iO});1
                
                %Define the beam loads & direct stress
                bl = i_defineBeamLoadResp(OptFEM(iO) , optPartCode{iO});
                if Options.ExitAfterSensitivity
                    %Make dummy constraints for the beam load DRESP2
                    BeamLoadResp   = bl(:);
                    BeamLoadConstr = arrayfun(@(~) awi.fe.opt.DesignConstraint, 1 : numel(BeamLoadResp));
                    set(BeamLoadConstr, {'DesignResponse'}, num2cell(BeamLoadResp));
                    addFEData(OptFEM(iO), bl, BeamLoadConstr);
                    continue
                end
                ds = i_defineBeamDirStress(OptFEM(iO), optPartCode{iO});
                addFEData(OptFEM(iO), bl, ds);
                
                %Define responses relating to the strength constraint set
                if Options.StrengthConstr
                    %Make the DRESP2/DCONSTR objects
                    ba  = i_defineBoxInternalArea(OptFEM(iO), optPartCode{iO});
                    ss  = i_defineBeamShrStress(OptFEM(iO)  , optPartCode{iO});
                    vm  = i_defineBeamVMStress(OptFEM(iO)   , optPartCode{iO});
                    rf  = i_defineBeamRF(OptFEM(iO), optPartCode{iO});
                    rfc = i_defineBeamRFConstr(Options);
                    %Add to the optimisation model
                    addFEData(OptFEM(iO), ba, ss, vm, rf, rfc);
                end
                
                %Define resposnes relating to stability of stringers/panels
                if Options.LocalStabConstr
                    %Make the DRESP2/DCONSTR objects
                    [ts, sa, tstr, zstr, si] = i_defineSkinAndStringerGeometry(OptFEM(iO), optPartCode{iO});
                    pb               = i_definePanelBuckling(optPartCode{iO});
                    sb               = i_defineStringerBuckling(optPartCode{iO});
                    [pb_rf , sb_rf]  = i_defineBucklingRF(optPartCode{iO});
                    [pb_rfc, sb_rfc] = i_definePanelBuckConstr;
                    %Add to the optimisation model
                    addFEData(OptFEM(iO), ts, sa, tstr, zstr, si, pb, sb, ...
                        pb_rf, sb_rf, pb_rfc, sb_rfc);
                end
                
            end
            
            if Options.EulerBuckConstr && ~Options.ExitAfterSensitivity
                %   - TODO: This needs to be generalised or replaced by the
                %   linear buckling analysis.
                euler_buckling(OptFEM, optPartCode, Options.GlobalBuckRF);
            end
            
            %Define the constraint set for strength & stability constraints
            dcon = [OptFEM.DesConstr];
            if ~isempty(dcon)
                DesConstr = awi.fe.opt.DesignConstraintSet;
                DesConstr.DesignConstraints = dcon;
                addFEData(FEM(1), DesConstr);
                ConstraintData.Strength = DesConstr;
            end
            
            function stbdInd = returnStarboardIndex(optPartNames, ind)
                
                %Determine which part we are dealing with (Wing/Strut/Jury/etc.)
                part = optPartNames{ind};
                %Find the part token
                temp = strsplit(part, 'Port');
                temp(contains(temp, 'Port')) = [];
                part = strtrim(strjoin(temp));
                %Find the 'Starboard' element and grab handle to those design vars
                stbdInd = find(contains(optPartNames, ['Starboard ', part]));
                
            end
            
            function [hLab, wLab] = i_grabBoxHeightAndWidth(OptFEM, optPartCode)
                %i_grabBoxHeightAndWidth Returns the labels for the box
                %height and width from the design table in 'OptFEM'.
                
                hLab = [];
                wLab = [];
                
                hTok = sprintf('%s_H', optPartCode);
                wTok = sprintf('%s_W', optPartCode);
                
                %Find the design table containing the data
                dt   = [OptFEM.DesTable];
                dt   = dt(arrayfun(@(dt) any(contains(dt.Prop, hTok)), dt));
                
                if isempty(dt)
                    return
                end
                
                %Labels of box weight/width in the Design Table
                wLab = dt.Prop(contains(dt.Prop, wTok));
                hLab = dt.Prop(contains(dt.Prop, hTok));
                
            end
            
            function matProps = i_grabMaterialProps(OptFEM, optPartCode)
                %i_grabMaterialProps Returns the labels for the material
                %properties from the design table in 'OptFEM'.
                
                matProps = [];
                
                matTok = sprintf('RHO_%s', optPartCode);
                
                %Find the design table containing the data
                dt   = [OptFEM.DesTable];
                dt   = dt(arrayfun(@(dt) any(contains(dt.Prop, matTok)), dt));
                
                if isempty(dt)
                    return
                end
                
                %Return all properties
                matProps = dt.Prop;
                
            end
            
            function ba = i_defineBoxInternalArea(OptFEM, optPartCode)
                %i_defineBoxInternalArea Defines the DRESP2 entry that
                %calculates the internal area of the box.
                
                ba = arrayfun(@(~) awi.fe.opt.DesignResponseEqn, 1 : nBP + 1);
                dv = OptFEM.DesVar;
                
                [hLab, wLab] = i_grabBoxHeightAndWidth(OptFEM, optPartCode);
                dims  = [wLab ; hLab];
                
                set(ba, 'DesignEquation', BoxArea);
                lab = arrayfun(@(ii) [optPartCode, 'A', num2str(ii)], 1 : nBP + 1, 'Unif', false);
                set(ba, {'LABEL'}, lab');
                set(ba, {'DesignVariables'}, num2cell(reshape(dv, [nOptVar, nBP + 1])', 2));
                set(ba, {'DesignConstants'}, num2cell(dims', 2));
                
            end
            
            function bl = i_defineBeamLoadResp(OptFEM, optPartCode)
                %i_defineBeamLoadResp Defines the DRESP1 entry for the 6
                %beam loads at ends A and ends B (final beam element only).
                
                nLoad = 6;
                nBP   = numel(OptFEM.BeamProps);
                nBR   = (nBP + 1) * nLoad;
                
                %Make the 'DesignResponse' object
                bl    = arrayfun(@(~) awi.fe.opt.DesignResponse, 1 : nBR);
                n   = numel(bl) / nLoad;
                bl  = reshape(bl, [nLoad, n]);
                
                set(bl, 'ResponseType', 'FORCE');
                set(bl, 'PropertyType', 'PBEAM');
                
                %Define the label
                %                 lab = strcat(repmat({optPartCode}, [6, n]), repmat({'M1' ; 'M2' ; 'F1' ; 'F2' ; 'AF' ; 'TQ'}, [1, n]));
                lab = strcat(repmat({optPartCode}, [6, n]), repmat({'Mi' ; 'Mo' ; 'Fi' ; 'Fo' ; 'AF' ; 'TQ'}, [1, n]));
                lab = strcat(lab, strcat('_', repmat(strsplit(num2str(1 : n)), [6, 1])));
                
                %Define the attributes
                itemCode = [4, 5, 6, 7, 8, 9];
                attr = cell(size(bl));
                for iCode = 1 : numel(itemCode)
                    %Assign item code for end-A
                    for iElem = 1 : nBP
                        attr{iCode, iElem} = {itemCode(iCode), [], OptFEM.BeamProps(iElem).ID};
                    end
                    %Assign item code for end-B
                    %    - Need to add (K-1) * 9 to the item code value, for end-B K=11
                    attr{iCode, end} = {itemCode(iCode) + 90, [], OptFEM.BeamProps(iElem).ID};
                end
                set(bl, {'ResponseAttributes'}, attr(:));
                set(bl, {'LABEL'}             , lab(:));
                
            end
            
            function ds = i_defineBeamDirStress(OptFEM, optPartCode)
                %i_defineBeamDirStress Defines the DRESP2 entries that
                %calculates the direct stress in the box.
                
                pr = OptFEM.DesPropRelEqn;
                
                ds = arrayfun(@(~) awi.fe.opt.DesignResponseEqn, 1 : nBP + 1);
                
                %Assign label
                lab = arrayfun(@(ii) [optPartCode, 'DS', num2str(ii)], 1 : nBP + 1, 'Unif', false);
                set(ds, {'LABEL'}, lab');
                
                %Assign DEQATN
                set(ds, 'DesignEquation', BoxStressEqn(1));
                
                %Assign DESVAR responses
                dv = OptFEM.DesVar;
                dv = reshape(dv, [nOptVar, nBP + 1]);
                set(ds, {'DesignVariables'}, num2cell(dv', 2));
                
                %Assign DTABLE properties
                arrayfun(@(i) set(ds(:, i), 'DesignConstants', dims(:, i)'), 1 : nBP + 1);
                
                %Assign DRESP1
                arrayfun(@(i) set(ds(:, i), 'DesignResponses', bl([1, 2, 5], i)'), 1 : nBP + 1);
                
                %Assign DVPREL2
                dvpr2   = reshape(pr, [nProp, nBP]);
                prpCode = {dvpr2(:, 1).PropNameOrField};
                idx     = ismember(prpCode, {'A(A)', 'I1(A)', 'I2(A)'});
                dvpr2A  = dvpr2(idx, :);
                dvpr2B  = dvpr2(ismember(prpCode, {'A(B)', 'I1(B)', 'I2(B)'}), end);
                dvpr2   = [dvpr2A, dvpr2B];
                arrayfun(@(i) set(ds(:, i), 'PropRelationEqn', dvpr2(:, i)'), 1 : nBP + 1);
                
            end
            
            function ss = i_defineBeamShrStress(OptFEM, optPartCode)
                %i_defineBeamShrStress Defines the DRESP2 entries that
                %calculates the shear stress in the box.
                %
                %   - Two shear stress values per section. One for the
                %   cover and one for the web.
                %   - Design variables are ordered by skin thickness then
                %   cover thickness
                
                dv   = OptFEM.DesVar;
                
                ss = arrayfun(@(~) awi.fe.opt.DesignResponseEqn, 1 : (nBP + 1) * nOptVar);
                ss = reshape(ss, [nOptVar, nBP + 1]);
                
                %Assign label
                lab = arrayfun(@(ii) { ...
                    [optPartCode, 'SSC', num2str(ii)] ; ...
                    [optPartCode, 'SSS', num2str(ii)]}, 1 : nBP + 1, 'Unif', false);
                lab = horzcat(lab{:});
                set(ss, {'LABEL'}, lab(:));
                
                %Assign DEQATN
                set(ss, 'DesignEquation', BoxStressEqn(2));
                
                %Assign DESVAR
                set(ss, {'DesignVariables'}, num2cell(dv(:)));
                
                %Assign DTABLE properties
                arrayfun(@(i) set(ss(:, i), 'DesignConstants', [dims(:, 1)', {'KLOAD'}]), 1 : nBP + 1);
                
                %Assign DRESP1
                arrayfun(@(i) set(ss(1, i), 'DesignResponses', bl([3, 4, 6], i)'), 1 : nBP + 1); % Cover
                arrayfun(@(i) set(ss(2, i), 'DesignResponses', bl([4, 3, 6], i)'), 1 : nBP + 1); % Spar
                
                %Assign DRESP2
                temp = [ba ; ba];
                set(ss, {'DesignResponseEqns'}, num2cell(temp(:)));
                
            end
            
            function vm = i_defineBeamVMStress(~, optPartCode)
                %i_defineBeamVMStress Defines the DRESP2 entries that
                %calculate the Von Mises stress in the box.
                %
                %   - Two Von Mises stress values per section. One for the
                %   cover and one for the web.
                
                %Design Responses - Von Mises Stress
                vm = arrayfun(@(~) awi.fe.opt.DesignResponseEqn, 1 : (nBP + 1) * nOptVar);
                vm  = reshape(vm, [nOptVar, nBP + 1]);
                
                %Define the label
                lab = arrayfun(@(ii) ...
                    {[optPartCode, 'VMC', num2str(ii)] ; [optPartCode, 'VMS', num2str(ii)]}, ...
                    1 : nBP + 1, 'Unif', false);
                lab = horzcat(lab{:});
                set(vm, {'LABEL'}, lab(:));
                
                %Assign DRESP2
                set(vm, 'DesignEquation', VMEqn);
                set(vm(1, :), {'DesignResponseEqns'}, num2cell([ds ; ss(1, :)]', 2));   %Cover
                set(vm(2, :), {'DesignResponseEqns'}, num2cell([ds ; ss(2, :)]', 2));   %Spar
                
            end
            
            function rf = i_defineBeamRF(OptFEM, optPartCode)
                %i_defineBeamRF Defines the DRESP2 entries that calculate
                %the Reserve Factors in the box.
                %
                %   - Two RF values per section. One for the cover and one
                %   for the web.
                
                rf = arrayfun(@(~) awi.fe.opt.DesignResponseEqn, 1 : (nBP + 1) * nOptVar);
                rf  = reshape(rf, [nOptVar, nBP + 1]);
                
                %Define the label
                lab = arrayfun(@(ii) { ...
                    [optPartCode, 'CRF', num2str(ii)] ; ...
                    [optPartCode, 'SRF', num2str(ii)]}, 1 : nBP + 1, 'Unif', false);
                lab = horzcat(lab{:});
                set(rf, {'LABEL'}, lab(:));
                
                %Define the DRESP2
                set(rf, {'DesignResponseEqns'}, num2cell(vm(:)));
                
                %Design Constants
                idx = arrayfun(@(dt) any(contains(dt.Prop, 'STE')), OptFEM.DesTable);
                dt = OptFEM.DesTable(idx);
                idx = contains(dt.Prop, 'STE');
                set(rf, 'DesignConstants'  , dt.Prop(idx));
                
                %Design Equation
                set(rf, 'DesignEquation', RFEqn);
                
            end
            
            function rfc = i_defineBeamRFConstr(Options)
                %i_defineBeamRFConstr Defines the DCONSTR entry for the
                %cover and spar web strength constraints.
                
                rfc = arrayfun(@(~) awi.fe.opt.DesignConstraint, 1 : (nBP + 1) * nOptVar);
                set(rfc, {'DesignResponse'}, num2cell(rf(:)));
                set(rfc, 'LowerBound', Options.StrengthRF);
                set(rfc, 'UpperBound', Options.MaxRF);
                
            end
            
            function [ts, sa, tstr, zstr, si]  = i_defineSkinAndStringerGeometry(OptFEM, optPartCode)
                %i_defineSkinAndStringerGeometry Defines the DRESP2 that
                %calculates the skin thickness (ts), stringer area (sa),
                %stringer neutral axis (zstr) and stringer second moment of
                %area (si).
                
                n = nBP + 1;
                
                %Grab design variables
                dv = OptFEM.DesVar;
                dv = reshape(dv, [nOptVar, n]);
                
                %Only want design variables related to the skin
                dv = dv(1, :);
                
                %Make the objects
                ts   = arrayfun(@(~) awi.fe.opt.DesignResponseEqn, 1 : n);
                sa   = arrayfun(@(~) awi.fe.opt.DesignResponseEqn, 1 : n);
                tstr = arrayfun(@(~) awi.fe.opt.DesignResponseEqn, 1 : n);
                zstr = arrayfun(@(~) awi.fe.opt.DesignResponseEqn, 1 : n);
                si   = arrayfun(@(~) awi.fe.opt.DesignResponseEqn, 1 : n);
                
                %Define label
                lab = arrayfun(@(i) sprintf('%sTSK%i', optPartCode, i), 1 : n, 'Unif', false);
                set(ts, {'LABEL'}, lab');
                lab = arrayfun(@(i) sprintf('%sAS%i', optPartCode, i), 1 : n, 'Unif', false);
                set(sa, {'LABEL'}, lab');
                lab = arrayfun(@(i) sprintf('%sTSt%i', optPartCode, i), 1 : n, 'Unif', false);
                set(tstr, {'LABEL'}, lab');
                lab = arrayfun(@(i) sprintf('%sZSt%i', optPartCode, i), 1 : n, 'Unif', false);
                set(zstr, {'LABEL'}, lab');
                lab = arrayfun(@(i) sprintf('%sIS%i', optPartCode, i), 1 : n, 'Unif', false);
                set(si, {'LABEL'}, lab');
                
                %Assign design equation
                set(ts  , 'DesignEquation', SkinThick);
                set(sa  , 'DesignEquation', StrArea);
                set(tstr, 'DesignEquation', StrThick);
                set(zstr, 'DesignEquation', StrZna);
                set(si  , 'DesignEquation', StrMomArea);
                
                %Assign DESVAR
                set(ts, {'DesignVariables'}, num2cell(dv'));
                set(sa, {'DesignVariables'}, num2cell(dv'));
                
                %Assign Design Constants
                set(ts  , 'DesignConstants', {'AS_BT'});
                set(sa  , 'DesignConstants', {'B_STR', 'AS_BT'});
                set([tstr, zstr, si], 'DesignConstants'  , {'W_STR'});
                
                %Assign DRESP2
                set(tstr, {'DesignResponseEqn'}, num2cell([sa ; ts]'  , 2));
                set(zstr, {'DesignResponseEqn'}, num2cell([ts ; tstr]', 2));
                set(si  , {'DesignResponseEqn'}, num2cell([ts ; tstr ; zstr]', 2));
                
            end
            
            function pb  = i_definePanelBuckling(optPartCode)
                %i_definePanelBuckling Defines the DRESP2 that calculates
                %the panel buckling stress.
                
                pb = arrayfun(@(~) awi.fe.opt.DesignResponseEqn, 1 : nBP + 1);
                
                %Define label
                lab = arrayfun(@(i) sprintf('%sPB%i', optPartCode, i), 1 : numel(pb), 'Unif', false);
                set(pb, {'LABEL'}, lab');
                
                %Assign design equation
                set(pb, 'DesignEquation', PanBuckEqn);
                
                %Assign design constants
                prps = [ ...
                    {'PI', 'KC', 'B_STR'}, ...
                    strcat({'E', 'NU'}, ['_', optPartCode])];
                set(pb, 'DesignConstants', prps);
                
                %Assign DRESP2
                set(pb, {'DesignResponseEqns'}, num2cell(ts)');
                
            end
            
            function sb  = i_defineStringerBuckling(optPartCode)
                %i_definePanelBuckling Defines the DRESP2 that calculates
                %the stringer buckling stress.
                
                sb = arrayfun(@(~) awi.fe.opt.DesignResponseEqn, 1 : nBP + 1);
                
                %Define label
                lab = arrayfun(@(i) sprintf('%sSB%i', optPartCode, i), 1 : numel(pb), 'Unif', false);
                set(sb, {'LABEL'}, lab');
                
                %Assign design equation
                set(sb, 'DesignEquation', StrBuckEqn);
                
                %Assign design constants
                prps = {'PI', ['E_', optPartCode], 'L_RIB'};
                set(sb, 'DesignConstants', prps);
                
                %Assign DRESP2
                set(sb, {'DesignResponseEqns'}, num2cell([sa ; si]', 2));
                
            end
            
            function [pb_rf, sb_rf] = i_defineBucklingRF(optPartCode)
                %i_defineBucklingRF Defines the DRESP2 entries that
                %calculate the Reserve Factors for panel and stringer
                %buckling.
                
                n = nBP + 1;
                
                pb_rf = arrayfun(@(~) awi.fe.opt.DesignResponseEqn, 1 : n);
                sb_rf = arrayfun(@(~) awi.fe.opt.DesignResponseEqn, 1 : n);
                
                %Add the label
                lab = arrayfun(@(i) sprintf('%sPRF%i', optPartCode, i), 1 : n, 'Unif', false);
                set(pb_rf, {'LABEL'}, lab');
                lab = arrayfun(@(i) sprintf('%sBRF%i', optPartCode, i), 1 : n, 'Unif', false);
                set(sb_rf, {'LABEL'}, lab');
                
                %Set design equation
                set([pb_rf, sb_rf], 'DesignEquation', RFEqn);
                
                %Define the DRESP2
                set(pb_rf, {'DesignResponseEqns'}, num2cell([pb ; ds(1, :)]', 2));
                set(sb_rf, {'DesignResponseEqns'}, num2cell([sb ; ds(1, :)]', 2));
                
            end
            
            function [pb_rfc, sb_rfc] = i_definePanelBuckConstr
                %i_definePanelBuckConstr Defines the DCONSTR entry for the
                %panel and stringer buckling constaints.
                
                pb_rfc = arrayfun(@(~) awi.fe.opt.DesignConstraint, 1 : (nBP + 1));
                sb_rfc = arrayfun(@(~) awi.fe.opt.DesignConstraint, 1 : (nBP + 1));
                set(pb_rfc, {'DesignResponse'}, num2cell(pb_rf(:)));
                set(sb_rfc, {'DesignResponse'}, num2cell(sb_rf(:)));
                set([pb_rfc, sb_rfc], 'LowerBound', Options.StabilityRF);
                set([pb_rfc, sb_rfc], 'UpperBound', Options.MaxRF);
                
            end
            
            %% Define Flutter responses & constraints
            
            if Options.FlutterConstr
                %Make a design response which will refer to the flutter damping
                FlutterResp = awi.fe.opt.DesignResponse;
                FlutterResp.LABEL        = 'FLTRDAMP';
                FlutterResp.ResponseType = 'FLUTTER';
                FlutterResp.PropertyType = 'PKNL';
                FlutterResp.ResponseAttributes = {[], [], [FlutterData.SET1.SID, FlutterData.FLFACT.SID]};
                %Make a design constraint that references this response
                FConstr = awi.fe.opt.DesignConstraint;
                FConstr.DesignResponse = FlutterResp;
                FConstr.LowerBound     = -2;
                FConstr.UpperBound     = 1e-10;
                %Add the data to the optimisation model
                addFEData(FEM(1), FlutterResp, FConstr);
                %Update the constraint data
                ConstraintData.Flutter = FConstr;
            end
            
            %% Define global buckling responses & constraints
            
            if Options.GlobalBuckConstr
                %Make a design response which will refer to the buckling eigenvalue
                GlobalBuck = awi.fe.opt.DesignResponse;
                GlobalBuck.LABEL              = 'BUCKEIG';
                GlobalBuck.ResponseType       = 'LAMA';
                GlobalBuck.ResponseAttributes = {1, [], []};    %Only interested in the first mode
                %Make a design constraint that references this response
                GBConstr = awi.fe.opt.DesignConstraint;
                GBConstr.DesignResponse = GlobalBuck;
                GBConstr.LowerBound     = Options.GlobalBuckRF;
                GBConstr.UpperBound     = Options.MaxRF;
                %Add the data to the optimisation model
                addFEData(FEM(1), GlobalBuck, GBConstr);
                %Update the constraint data
                ConstraintData.GlobalBuckling = GBConstr;
            end
            
            %% Export FEM & OptM to files
            
            %Update the ID numbers
            assignIDnumbers([allFEM, MassModel]);
            
            %Export data?
            if p.Results.WriteDataFiles
                
                %Write the model bulk data
                [~, partNames] = export(FEM, outdir, ...
                    'CollectChildModels'   , true  , ...
                    'FullyQualifiedPath'   , false , ...
                    'WriteHeaderFile'      , false , ...
                    'WriteOptimisationData', false , ...
                    'WriteDesignModel'     , true);
                
                %Write the model mass cases (seperate files)
                massNames = cell(size(MassModel));
                for iM = 1 : numel(massNames)
                    [~, massNames(iM)] = export(MassModel(iM), outdir, ...
                        'FullyQualifiedPath'   , false , ...
                        'WriteHeaderFile'      , false );
                end
                
                %Write the optimisation model
                %   - Split the data by :
                %       * Design Variables
                desOptFiles{1} = 'DesVar.bdf';
                fileName       = fullfile(outdir, desOptFiles{1});
                fid            = fopen(fileName, 'w');
                writeToFile([allFEM.DesVar], fid);
                fclose(fid);
                %       * Design Relations
                desOptFiles{2} = 'DesRel.bdf';
                fileName       = fullfile(outdir, desOptFiles{2});
                fid            = fopen(fileName, 'w');
                writeToFile([allFEM.DesPropRel]   , fid);
                writeToFile([allFEM.DesPropRelEqn], fid);
                fclose(fid);
                %       * Design Responses
                desOptFiles{3} = 'DesResp.bdf';
                fileName       = fullfile(outdir, desOptFiles{3});
                fid            = fopen(fileName, 'w');
                writeToFile([allFEM.DesResp]   , fid);
                writeToFile([allFEM.DesRespEqn], fid);
                fclose(fid);
                %       * Design Constraints
                desOptFiles{4} = 'DesConstr.bdf';
                fileName       = fullfile(outdir, desOptFiles{4});
                fid            = fopen(fileName, 'w');
                writeToFile([allFEM.DesConstr], fid);
                writeToFile([allFEM.DesConstrSet], fid);
                fclose(fid);
                %       * Design Equations and Tables
                desOptFiles{5} = 'DesEqn.bdf';
                fileName       = fullfile(outdir, desOptFiles{5});
                fid            = fopen(fileName, 'w');
                writeToFile([allFEM.DesTable], fid);
                writeToFile([allFEM.DesEqn]  , fid);
                fclose(fid);
                
            else
                
                partNames   = [];
                massNames   = [];
                desOptFiles = [];
                
            end
            
            %% Save the FEM objects & inputs to .mat files for post-processing
            
            RunData.Options  = Options;
            RunData.PartName = optPartName;
            RunData.PartCode = optPartCode;
            RunData.VarNames = optVarNames;
            RunData.VarLabel = optLabel;    %#ok<STRNU>
            
            save(fullfile(outdir, 'OptData.mat'), 'FEM', 'MassModel');
            save(fullfile(outdir, 'RunData.mat'), 'RunData');
            
            %% Write the SOL200 header file

            if p.Results.WriteHeaderFile && p.Results.WriteDataFiles
                %Files to be included in the analysis
                includeFiles = [partNames, desOptFiles, massNames];
                %Write the header
                optFile = writeSol200HeaderFile(obj, Aircraft        , ...
                    LoadCases, Options, outdir, TrimData, StaticData , ...
                    BuckData, FlutterData, WeightResp, ConstraintData, ...
                    includeFiles, 'MassModels', MassSets             , ...
                    'SensitivityAnalysis', Options.ExitAfterSensitivity);
            else
                optFile = [];
            end
            
            if nargout > 0
                varargout{1} = FEM;
            end
            if nargout > 1
                varargout{2} = OptFEM;
            end
            if nargout > 2
                varargout{3} = partNames;
            end
            if nargout > 3
                varargout{4} = desOptFiles;
            end
            if nargout > 4
                varargout{5} = massNames;
            end
            
        end
        function optFile = writeSol200HeaderFile(obj, Aircraft, LoadCases, Options, outdir, TrimData, StaticData, BuckData, FlutterData, WeightResp, ConstraintData, includeFiles, varargin)
            %writeSol200HeaderFile Writes the main header file for running
            %the SOL 200 design optimisation.
            
            p = inputParser;
            addParameter(p, 'MassModels'         , []   );
            addParameter(p, 'WriteAnalysisData'  , true , @(x)validateattributes(x, {'logical'}, {'scalar'}));
            addParameter(p, 'FullyQualifiedPath' , false, @(x)validateattributes(x, {'logical'}, {'scalar'}));
            addParameter(p, 'SensitivityAnalysis', false, @(x)validateattributes(x, {'logical'}, {'scalar'}));
            parse(p, varargin{:});
            
            optFile = fullfile(outdir, 'SOL200_header.bdf');
            fid = fopen(optFile, 'w');
            
            %What load cases have we got?
            LCType   = get(LoadCases, {'LoadCaseType'});
            idxT     = ismember(LCType, {'Manoeuvre', 'Pratt Gust'});
            idxS     = ismember(LCType, {'Static'});
            bTrim    = any(idxT);
            bStatic  = any(idxS);
            bFlutter = [];
            
            %% File Management Statement (FMS)
            %             if Options.GlobalBuckConstr
            %                 awi.fe.FEBaseClass.writeHeading(fid, 'F I L E  M A N A G E M E N T  S T A T E M E N T');
            %                 %   - Define the various output files for trim loads etc.
            %                 fileNo = 53;
            %                 fileNo = fileNo : fileNo + numel(LoadCases) - 1;
            %                 trimForceFiles   = arrayfun(@(i) sprintf('trimLoadsSC%i.inc', i), 1 : numel(fileNo), 'Unif', false);
            %                 trimForceCommand = arrayfun(@(i) sprintf(['ASSIGN USERFILE=', ...
            %                     '''%s'',STATUS=UNKNOWN,FORMATTED,UNIT=%i'], ...
            %                     trimForceFiles{i}, fileNo(i)), 1 : numel(fileNo), 'Unif', false);
            %                 fprintf(fid, '%-s\r\n', trimForceCommand{:});
            %                 %   - Make dummy files so that Nastran run does not fail when processing
            %                 %   the header file.
            %                 for jj = 1 : numel(trimForceFiles)
            %                     temp = fopen(trimForceFiles{jj}, 'w');
            %                     fclose(temp);
            %                 end
            %                 clear temp jj
            %             else
            %                 trimForceFiles = [];
            %                 fileNo = [];
            %             end
            
            trimForceFiles = [];
            fileNo = [];
            
            %% Executive Control
            awi.fe.FEBaseClass.writeHeading(fid, 'E X E C U T I V E  C O N T R O L');
            fprintf(fid, 'SOL %i\r\n', 200);
            fprintf(fid, 'ECHOOFF         $ SUPPRESSES THE ECHO OF EXECUTIVE CONTROL\r\n');
            fprintf(fid, 'CEND\r\n');
            
            %% Case Control
            obj.defaultCaseControl(fid, false); %Printing the No. lines seems to cause an error...
            %   - SOL 200
            fprintf(fid, '$O B J E C T I V E   F U N C T I O N\r\n$ - Minimise weight\r\n');
            fprintf(fid, 'DESOBJ(MIN) = %i\r\n', WeightResp.ID);
            if p.Results.SensitivityAnalysis    %Sensitivity analysis?
                fprintf(fid, ['$ * * Analysis will exit after sensitivities ', ...
                    'have been calculated * * \nDSAPRT(', ...
                    'FORMATTED,NOEXPORT,START=1,END=SENS) = ALL\r\n']);
            else
                fprintf(fid, 'DSAPRT(FORMATTED,NOEXPORT,START=1,BY=1,END=LAST) = ALL\r\n'); %Request sensitivities
            end
            pad = blanks(4);
            %   - SOL 144
            sc0 = 1;
            if bTrim
                obj.printTrimSubcases(fid, sc0, TrimData, ConstraintData.Strength, pad, fileNo);
            end
            %   - SOL 101
            if bStatic
                sc0 = sc0 + numel(TrimData.Trim);
                obj.printStaticSubcases(fid, sc0, StaticData, ConstraintData.Strength, pad);
            end
            %   - SOL 105
            if Options.GlobalBuckConstr
                obj.printBucklingSubcases(fid, TrimData, pad, ConstraintData.GlobalBuckling, GBConstr);
            end
            %   - SOL 145
            if Options.FlutterConstr
                sc0 = sc0 + 1;
                obj.printFlutterSubcases(fid, FlutterData, ConstraintData.Flutter, sc0, pad);
            end
            
            %% Bulk Data
            obj.defaultBulkStatement(fid);
            %   - Optimisation parameters
            obj.printDesOptParam(fid, Options);
            if ~Options.ScreenConstr            %Constraint screening?
                if Options.StrengthConstr || Options.LocalStabConstr
                    fprintf(fid, '%-8s%-8s%-8s%-8i\r\n', 'DSCREEN', 'EQUA' , '-1.e+12', 5e4);
                    fprintf(fid, '%-8s%-8s%-8s%-8i\r\n', 'DSCREEN', 'FORCE', '-1.e+12', 5e4);
                end
                if Options.FlutterConstr
                    fprintf(fid, '%-8s%-8s%-8s%-8i\r\n', 'DSCREEN', 'FLUTTER', '-1.e+12', 5e4);
                end
            end
            if p.Results.WriteAnalysisData      %Write bulk data for analysis or in external file?
                %fprintf(fid, 'DOPTPRM,P1,1,P2,14\r\n'); %Optimisation diagnostics
                %   - Trim Data
                if bTrim
                    TrimLoadCases = LoadCases(idxT);
                    %Open the file
                    trimFile = fullfile(outdir, 'TrimData.bdf');
                    tfid     = fopen(trimFile, 'w');
                    awi.fe.FEBaseClass.writeColumnDelimiter(tfid, 'normal');
                    %spc
                    fprintf(tfid, '%-8s%-8i%-8i%-8i\r\n', 'SPC1', TrimData.SPC_id, TrimData.SPCdof, TrimData.RefGrid.ID);
                    %suport
                    fprintf(tfid, '%-8s%-8i%-8i\r\n', 'SUPORT', TrimData.RefGrid.ID, TrimData.SUPdof);
                    %aeros
                    obj.writeSteadyAeroEntry(tfid, Aircraft, 'SymXZ', Options.ModelXZSymmetry);
                    %trim
                    obj.writeTrimEntries(tfid, TrimData.Trim, TrimLoadCases);
                    %aestat
                    data = [num2cell(TrimData.AESTAT_id) ; TrimData.AESTAT];
                    str  = sprintf('AESTAT,%i,%s-', data{:});
                    data = strsplit(str, '-');
                    fprintf(tfid, '%s\r\n', data{:});
                    fprintf(tfid, 'PARAM,AUNITS,0.1019716\r\n');
                    fclose(tfid);
                else
                    trimFile = [];
                end
                %   - Static Data
                if bStatic
                    StaticLoadCases = LoadCases(idxS);
                    %Open the file
                    staticFile = fullfile(outdir, 'StaticData.bdf');
                    sfid       = fopen(staticFile, 'w');
                    awi.fe.FEBaseClass.writeColumnDelimiter(sfid, 'normal');
                    %spc
                    fprintf(sfid, '%-8s%-8i%-8i%-8i\r\n', 'SPC1', ...
                        StaticData.SPC_id, StaticData.SPCdof, StaticData.RefGrid.ID);
                    %grav
                    fprintf(sfid, '%-8s%-8i%-8i%-8.5f%-8.3f%-8.3f%-8.3f\r\n', ...
                        'GRAV', StaticData.GRAV_id, 0, 9.80665, 0, 0, -1);
                    %load
                    for i = 1 : numel(StaticLoadCases)
                        fprintf(sfid, '%-8s%-8i%-8.3f%-8.3f%-8i\r\n', 'LOAD', ...
                            StaticData.LOAD_id(i), 1, StaticLoadCases(i).LoadFactor, ...
                            StaticData.GRAV_id);
                    end
                    fclose(sfid);
                else
                    
                    staticFile = [];
                end
                %   - Buckling Data
                if Options.GlobalBuckConstr
                    buckFile = fullfile(outdir, 'BucklingData.bdf');
                    bfid = fopen(buckFile, 'w');
                    %spc
                    fprintf(bfid, '%-8s%-8i%-8i%-8i\r\n', 'SPC1', BuckData.SPC_id, BuckData.SPCdof, BuckData.RefGrid.ID);
                    %method
                    fprintf(bfid, '%-8s%-8i%-#8.3g%-8s%-8i\r\n', 'EIGRL', BuckData.EIGRL_id, 0, blanks(8), BuckData.NModes);
                    buckFile = {buckFile};
                    fclose(bfid);
                else
                    buckFile = [];
                end
                %   - Flutter Data
                if Options.FlutterConstr
                    flutterFile = fullfile(outdir, 'FlutterData.bdf');
                    ffid = fopen(flutterFile, 'w');
                    awi.fe.FEBaseClass.writeColumnDelimiter(ffid, 'normal');
                    %spc
                    fprintf(ffid, '%-8s%-8i%-8i%-8i\r\n', 'SPC1', FlutterData.SPC.SID, FlutterData.SPC.C, FlutterData.SPC.G);
                    %method
                    fprintf(ffid, '%-8s%-8i%-#8.3g%-8s%-8i\r\n', 'EIGRL', FlutterData.EIGRL.SID, 0, blanks(8), FlutterData.EIGRL.ND);
                    %aero
                    obj.writeUnsteadyAeroEntry(ffid, Aircraft, LoadCases, 'SymXZ', Options.ModelXZSymmetry);
                    %flutter
                    fprintf(ffid, '%-8s%-8i%-8s%-8i%-8i%-8i%-8s%-8s%-8s\r\n', 'FLUTTER', ...
                        FlutterData.FLUTTER.SID  , FlutterData.FLUTTER.METHOD, ...
                        FlutterData.FLUTTER.DENS , FlutterData.FLUTTER.MACH  , ...
                        FlutterData.FLUTTER.RFREQ, 'L', blanks(8), '1.E-4');
                    %flfact
                    for i = 1 : numel(FlutterData.FLFACT)
                        n = numel(FlutterData.FLFACT(i).F);
                        fmt = ['%-8s%-8i', repmat('%-8s', [1, 8]), '\r\n'];
                        %         fmt = ['%-8s%-8i', repmat('%-8.3f', [1, n]), repmat('%-8s', [1, 7-n]), '%-8s\r\n'];
                        F_  = cellstr(num2str(FlutterData.FLFACT(i).F, '%-8.3f'))';
                        dat = [{'FLFACT', FlutterData.FLFACT(i).SID}, F_, repmat({blanks(8)}, [1, 7-n]), {FlutterData.FLFACT(i).LAB}];
                        fprintf(ffid, fmt, dat{:});
                    end
                    %mkaero
                    m = num2cell(FlutterData.MKAERO.M);
                    k = num2cell(FlutterData.MKAERO.K);
                    if numel(k) > 8 || numel(m) > 8
                        error('Update code for MKAERO entry with more than 8 K/M terms');
                    end
                    fprintf(ffid, ['%-8s', repmat('%#-8.3g', [1, numel(m)]), '\r\n'], 'MKAERO1', m{:});
                    fprintf(ffid, ['%-8s', repmat('%#-8.3g', [1, numel(k)]), '\r\n'], blanks(8), k{:});
                    %set1
                    fprintf(ffid, '%-8s%-8i%-8i%-8s%-8i\r\n', 'SET1', FlutterData.SET1.SID, ...
                        FlutterData.SET1.ID(1), 'THRU', FlutterData.SET1.ID(end));
                    %tabdmp1
                    dat = [FlutterData.TABDMP1.x ; FlutterData.TABDMP1.y];
                    str = cellstr(num2str(dat(:), '%#-8.2g'));
                    str = [str ; {'ENDT'}];
                    if numel(str) > 8
                        error('Update code for writing structural damping table');
                    end
                    fprintf(ffid, ['%-8s%-8i%-8s\r\n%-8s', repmat('%-8s', [1, numel(str)]), '\r\n'], ...
                        'TABDMP1', FlutterData.TABDMP1.TID, 'CRIT', blanks(8), str{:});
                    fclose(ffid);
                    flutterFile = {flutterFile};
                else
                    flutterFile = [];
                end
            else
                trimFile    = [];
                staticFile  = [];
                buckFile    = [];
                flutterFile = [];
            end
            
            %Mass Sets
            if ~isempty(p.Results.MassModels)
                MassModel = p.Results.MassModels;
                if isa(MassModel, 'awi.fe.MassModel')
                    writeToFile([MassModel.MassGroups], fid);
                elseif isa(MassModel, 'awi.fe.MassSet')
                    writeToFile(MassModel, fid);
                end
            end
            
            %Additional bulk data files
            includeFiles = [trimFile, staticFile, buckFile, flutterFile, trimForceFiles, includeFiles];
            if ~iscell(includeFiles)
                includeFiles = {includeFiles};
            end
            awi.fe.FEBaseClass.writeSubHeading(fid, 'I N C L U D E  F I L E S');
            awi.fe.FEBaseClass.writeIncludeStatement(fid, includeFiles, ...
                'FullyQualifiedPath', p.Results.FullyQualifiedPath);
            
            %End of file
            fprintf(fid, 'ENDDATA\r\n');
            
            %Close the file
            fclose(fid);
            
        end
        function bOpt = canOptimise(obj)
            %canOptimise Checks whether the user has provided enough
            %information to run a structural optimisation.
            
            bOpt = true;
            
        end
    end
    
    methods % Running Nastran via command line
        function runNastran(obj, datfile, varargin)
            %runNastran Runs the MSC.Nastran executable with the filename
            %specified by 'datfile'.
            %
            % Details:
            %   - The 'system' command is used to run Nastran.
            %   - Must change the current directory to the local directory
            %     of the .dat file otherwise the output files will be sent
            %     to the current working directory.
            %   - The 'system' command does not wait for the program to
            %     finish. Instead, we pause Matlab and wait for the
            %     analysis to end.
            %   - Whilst the .f06 is being checked for the end statement,
            %     it is also parsed for FATAL errors.
            
            p = inputParser;
            addParameter(p, 'MonitorProgress', 'Default', @(x)any(validatestring(x, {'Default', 'Nonlinear'})));
            parse(p, varargin{:});
            
            %Remember where the analysis is being invoked from
            this_dir = pwd;
            bFinish  = false;
            
            %Get the path to the MSC.Nastran executable
            nastran_exe = obj.NastranExe;
            if isempty(nastran_exe)
                return
            end
            
            %Set up file information
            [path, name, ext] = fileparts(datfile);
            
            %Run MSC.Nastran - TODO: Add 'scratch' options as input
            cd(path)
            cmd = strjoin({nastran_exe, [name ext], 'scratch=yes'}, ' ');
            system(cmd);
            
            %Define name of .f06 file
            f06name = [lower(name), '.f06'];
            
            %Monitor progress?
            switch p.Results.MonitorProgress
                case 'Default'
                    logFcn = @(f06Line, Data) i_checkForFatal(f06Line, Data);
                    Data = [];
                case 'Nonlinear'
                    logFcn = @(f06Line, Data) i_displayNonlinearConvergence(f06Line, Data);
                    hF     = figure('Name', 'Nonlinear Solution Sequency Progress');
                    hAx(1) = axes('Parent', hF, 'NextPlot', 'add');
                    Data   = struct( ...
                        'SolutionData'        , struct('NLIter', zeros(21, 1)), ...
                        'GraphicsHandles'     , [hF, hAx], ...
                        'ExtractIterationData', false);
            end
            
            
            %Check the .f06 for "FATAL ERROR" and "END OF JOB"
            while ~bFinish
                %Start a timer
                t0 = tic;
                t  = 0;
                %Pause Matlab until the .f06 is generated
                while isempty(dir(f06name))
                    if t > obj.F06MaxWait
                        break
                    end
                    pause(0.01);
                    t = toc(t0);
                end
                %Once the .f06 exists open the file and search for keywords
                if isempty(dir(f06name))
                    cd(this_dir);
                    error(['MSC.Nastran failed to start the analysis.', ...
                        '\n\n\t-%-12s: %s\n\t-%-12s: %s\n\n'], 'File', ...
                        [name, ext], 'Directory', path);
                else
                    fid = fopen(f06name, 'r');
                    while feof(fid) ~= 1
                        f06Line = fgets(fid);
                        [bFinish, bFatal, Data] = logFcn(f06Line, Data);
                        if bFatal
                            fclose(fid);
                            diary off
                            error(['* * * FATAL ERROR HAS OCCURED IN ', ...
                                'THE FILE %s * * *'], fullfile(path, f06name));
                        end
                        if bFinish
                            break
                        end
                    end
                    fclose(fid);
                end
            end
            
            
            function [bFinish, bFatal, Data] = i_checkForFatal(f06Line, Data)
                %i_checkForFatal Checks the f06Line for fatal errors and
                %the end of the analysis.
                
                bFinish = false;
                bFatal  = false;    
                
                if f06Line == -1
                    return
                end
                
                %Check for fatal error
                if contains(f06Line, 'FATAL MESSAGE')
                    %SOL 200 has a specific line that contains a
                    %FATAL error. Check this is not that.
                    if ~contains(f06Line, 'IF THE FLAG IS FATAL')
                        bFatal = true;
                    end
                end
                
                %Check for end of .f06 file
                if contains(f06Line, '* * * END OF JOB * * *')
                    bFinish = 1;
                end
                        
            end
            
            function [bFinish, bFatal, Data] = i_displayNonlinearConvergence(f06Line, Data)
                %i_displayNonlinearConvergence Plots the convergence
                %progress of a nonlinear analysis.
                
                bFinish = false;
                bFatal  = false;
                
                if f06Line == -1
                    return
                end
                
                %Pass it on
                [bFinish, bFatal, Data] = i_checkForFatal(f06Line, Data);
                if bFinish || bFatal
                    return
                end
                
                if contains(f06Line, 'N O N - L I N E A R   I T E R A T I O N   M O D U L E   O U T P U T')
                    Data.ExtractIterationData = true;
                end
                if Data.ExtractIterationData
                    iterData = sscanf(f06Line, '%%%f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f');
                    if ~isempty(iterData)
                        Data.SolutionData.NLIter = horzcat(Data.SolutionData.NLIter, iterData);
                    end
                end
                
                
            end
            
            %Return to the invoking directory
            cd(this_dir);
            
        end
    end
    
    methods (Access = private) %Helper functions
        function fName = parseFilename(obj, fName)
            %parseFilename Checks that the name of the file is a valid
            %string that is acceptable for MSC.Nastran analysis.
            %
            % Details:
            %   - Blank spaces will be removed
            %   - The list of invalid characters in provided in
            %     'obj.InvalidFileChar'.
            
            %Remove blanks
            fName = strrep(fName, ' ', '');
            
            %Loop through other forbidden tokens
            for ii = 1 : numel(obj.InvalidFileChar)
                fName = strrep(fName, obj.InvalidFileChar{ii}, '_');
            end
            
        end
        function ResultsObj = generateAWIResults(obj, NastranResults)
            %generateAWIResults Creates the 'awi.model.BeamResult' objects.
            
            ResultsObj = [];
            if nargin < 2 || isempty(NastranResults)
                return
            end
                        
            if isfield(NastranResults, 'h5data')
                ResultsObj = generateResultsFromH5(obj, NastranResults.h5data, ResultsObj);
            end
                        
            if isfield(NastranResults,'f06data')
                ResultsObj = generateResultsFromF06(obj, NastranResults.f06data, ResultsObj);
            end
            
        end
        function ResultsObj = generateResultsFromH5(obj, H5Data, ResultsObj)
            %generateResultsFromH5 Generates the ALENA results objects
            %using the results data in 'H5Data'.
            
            if ~isfield(H5Data.Raw.NASTRAN, 'INPUT')
                warning(['The model data for the Nastran model must '    , ...
                    'be contained in the HDF5 file in order to generate ', ...
                    'the results sets. Update MDLPRM,HDF5 and run the '  , ...
                    'analysis again.']);
            end
            resNames = fieldnames(H5Data.ResultSets);
            
            if any(strcmpi(resNames, 'displacement')) %Displacements
                DispResults = generateDisplacementResults(H5Data, obj.AnalysisModel);
                %ResultsObj = [ResultsObj, DispResults];
                ResultsObj.displacement = DispResults;
            end
            
            if any(strcmpi(resNames, 'stress')) %Stress
                StressResults = generateStressResults(H5Data, obj.AnalysisModel);
                ResultsObj.stress = StressResults;
            end
            
        end
        function ResultsObj = generateResultsFromF06(obj, F06Data, ResultsObj)
            %generateResultsFromH5 Generates the ALENA results objects
            %using the results data in 'F06Data'.
            
            % worth developping ? (charles)
            ResultsObj=[];
            
            error('Update code to extract results from F06 file');
            
        end
    end
    
    methods (Static) %Generic Nastran methods
        function NastranResults = extractNastranResults(datfile, varargin)
            %extractNastranResults Reads the results data from the
            %MSC.Nastran results files.
            %
            % Details:
            %   - Will read data from .h5 & .f06 files
            
            NastranResults = [];
            
            %Prompt user if no file is provided
            if nargin == 0
                datfile = getfile({ ...
                    '*.f06', 'F06 File' ; '*.h5', 'HDF5 File'}, ...
                    'Select a MSC.Nastran results file');
                if isempty(datfile), return; end %Escape route
                %Check the extension and modify inputs as appropriate
                [~, ~, ext] = fileparts(datfile);
                switch ext
                    case '.f06'
                        if ~any(ismember(varargin, 'ReadF06'))
                            varargin = [varargin, {'ReadF06', true}];
                        end
                end
            end
            
            %Parse inputs
            p = inputParser;
            addParameter(p, 'ReadHDF5', true, @(x)validateattributes(x, {'logical'},  {'scalar'}));
            addParameter(p, 'ReadF06' , false, @(x)validateattributes(x, {'logical'},  {'scalar'}));
            addParameter(p, 'Results' , {}   , @iscellstr);
            parse(p, varargin{:});
            
            %Do we need to do anything?
            if ~any([p.Results.ReadHDF5, p.Results.ReadF06])
                return
            end
            
            %Set up files to read
            [path, name, ~] = fileparts(datfile);
            
            %Output name is always lower case
            name_ = lower(name);
            
            %Read the HDF5 output
            if p.Results.ReadHDF5
                h5file = fullfile(path, [name_, '.h5']);
                [NastranResults.h5data.Raw, ~, NastranResults.h5data.ResultSets] = h5extract(h5file);
            end
            
            %Read the .F06 output
            if p.Results.ReadF06
                f06file = fullfile(path, [name_, '.f06']);
                NastranResults.f06data = f06extract(f06file);
            end
            
            
        end
        function defaultBulkStatement(fid, varargin)
            %defaultBulkStatement Writes the standard comments and
            %parameters for the model bulk data.
            %
            % Parameter Inputs:
            %   - 'RequestXBD' : Requests the generation of the MSC.Patran
            %   .xdb database.
            %   Defualt = true
            %
            % Nastran Parameters:
            %   'MDLPRM,HDF5,1'       - Generates a compressed HDF5 output
            %   'MDLPRM,OFFDEF,LROFF' - Enables "Large rotation offsets" for beam elements
            %   'PARAM,POST,0'        - Generates a MSC.Patran .xdb file
            %   'PARAM,GRDPNT,0'      - Generates the mass data
            %   'PARAM,BAILOUT,-1'    - Allows the analysis to run with
            %                           mechanisms
            
            p = inputParser;
            addParameter(p, 'RequestXDB', true, @(x)validateattributes(x, {'logical'}, {'scalar'}));
            parse(p, varargin{:});
            
            awi.fe.FEBaseClass.writeHeading(fid, 'B E G I N  B U L K');
            fprintf(fid, 'BEGIN BULK\r\n');
            %Request HDF5 & XDB output as standard
            fprintf(fid, 'MDLPRM,HDF5,1\r\n');
            %Enables "Large rotation offsets" for beam elements (WiA and Wib defines the offset directions)
            fprintf(fid, 'MDLPRM,OFFDEF,LROFF\r\n');
            if p.Results.RequestXDB
                fprintf(fid, 'PARAM,POST,0\r\n');
            end
            %Request Grid Point Weight Generator (GPWG)
            fprintf(fid, 'PARAM,GRDPNT,0\r\n');
            %Anticipate mechanisms because of the presence of joints
            fprintf(fid, 'PARAM,BAILOUT,-1\r\n');
            
        end
        function defaultCaseControl(fid, bLine)
            %defaultCaseControl Writes the standard case control options
            %for the 'Case Control Section'.
            %
            % Parameters:
            %   'LINE = 999999...' - Sets maximum value for printed page
            %                        length
            %   'ECHO = NONE'      - Turns off bulk data echo
            
            if nargin < 2
                bLine = true;
            end
            
            awi.fe.FEBaseClass.writeHeading(fid, 'C A S E  C O N T R O L');
            awi.fe.FEBaseClass.writeSubHeading(fid, 'O U T P U T  O P T I O N S');
            if bLine
                fprintf(fid, 'LINE = 99999999   $ SPECIFIES THE NUMBER OF LINES PER PRINTED PAGE\r\n');
            end
            fprintf(fid, 'ECHO = NONE       $ SUPPRESSES THE ECHO OF BULK DATA\r\n');
        end
        function outdir = makeDefaultAnalysisDirectory(analysisType, dirLoc)
            %makeDefaultAnalysisDirectory Creates a directory for the
            %Nastran analysis of type 'analysisType'.
            
            %Parse
            validateattributes(analysisType, {'char'}, {'row', 'nonempty'}, ...
                'makeDefaultAnalysisDirectory', 'analysis type');
            if nargin < 2
                dirLoc = pwd;
            end
            
            %Make the folder for the analysis
            dirName = [datestr(now, 30), '_', analysisType, '_analysis'];
            [bSuccess, ~, ~] = mkdir(dirLoc, dirName);
            assert(bSuccess == true, sprintf(['Unable to create ', ...
                'output directory for the %s analsysis.'], analysisType));
            outdir = fullfile(dirLoc, dirName);
            
        end
    end
    
    methods (Static) %Static Aeroelastic methods
        function checkAircraftTrimData(Aircraft)
            %checkAircraftTrimData Checks if the aircraft can actually be
            %trimmed.
            
            prp = {'RefSpan', 'RefChord', 'RefArea'};
            idx = cellfun(@isempty, get(Aircraft, prp));
            if any(idx)
                error(['The following reference quantities have not ', ...
                    'been defined for the Aircraft (%s):\n\t%s\n'], ...
                    Aircraft.Name, strjoin(prp(idx), ', '));
            end
            
        end
        function writeSteadyAeroEntry(fid, Aircraft, varargin)
            %writeSteadyAeroEntry Writes the equivalent MSC.Nastran AEROS
            %entry using the data in 'Aircraft'.
            
            p = inputParser;
            addParameter(p, 'SymXZ', false, @(x)validateattributes(x, {'logical'}, {'scalar'}));
            addParameter(p, 'SymXY', false, @(x)validateattributes(x, {'logical'}, {'scalar'}));
            parse(p, varargin{:});
            
            %Check for trim data
            awi.methods.Nastran.checkAircraftTrimData(Aircraft)
            
            %Write the data to the file
            fprintf(fid, '%-8s%-8i%-8i%-8.4f%-8.4f%-8.4f%-8i%-8i\r\n', ...
                'AEROS', 0, 0, Aircraft.RefChord, Aircraft.RefSpan   , ...
                Aircraft.RefArea, p.Results.SymXZ, p.Results.SymXY);
            
        end
        function TrimData = getTrimData(FEM, Aircraft, LoadCases, MassCases, ID0, RefGrid,varargin)
            %getTrimData Returns a MATLAB structure describing the basic
            %trim data for the load cases described in 'LoadCases'.
                        
            TrimData = struct();
            
            p = inputParser;
            p.addParameter('AoA', NaN);
            p.addParameter('ControlDOFs',@(x)contains(x,{'elev','ail', 'rud'}))
            parse(p, varargin{:});
            
            %Do we have any control surfaces?
            %   - The presence of control surfaces dictates whether we can
            %     trim for the pitch DOF.
            %   - TODO : Should check the normal vector for each control
            %   surface and make sure there is a component in the global-Z.
            % ControlSurf = Aircraft.ControlSurfaces;
            
            %Collapse the model hierachy
            controlsurfs = awi.methods.Nastran.getControlSurfs(FEM);

            %What aero controllers have been defined?
            controlsurfs = awi.methods.Nastran.getControlSurfs(FEM);
            I = ~p.Results.ControlDOFs(controlsurfs);
            staticControlSurf = controlsurfs(I);
            staticControlSurf = [staticControlSurf;num2cell(zeros(size(staticControlSurf)))];
            
            % get AESTAT entries
            if isnan(p.Results.AoA) 

                if isempty(any(~I))
                    dof = 3;
                    spc = 12456;
                    AESTAT = {'ANGLEA',[];'URDD3',LoadCases.LoadFactor};

                elseif any([LoadCases.AileronAngle,LoadCases.ElevatorAngle,LoadCases.RudderAngle])
                    dof = 456;
                    spc = 123;

                    ailer_surf=controlsurfs(cellfun(@(x)contains(x,'aile'),controlsurfs));
                    elev_surf=controlsurfs(cellfun(@(x)contains(x,'elev'),controlsurfs));
                    ruder_surf=controlsurfs(cellfun(@(x)contains(x,'rudder'),controlsurfs));

                    AESTAT = {'URDD4',0; 'URDD5',0; 'URDD6',0;...
                        ailer_surf{1}, LoadCases.AileronAngle;...
                        elev_surf{1}, LoadCases.ElevatorAngle;...
                        ruder_surf{1}, LoadCases.RudderAngle;...
                        'Roll', []; 'Pitch', []; 'Yaw', []};

                elseif ~any(LoadCases.AileronAngle)
                    dof = 35;
                    spc = 1246;
                    AESTAT = {'ANGLEA',[];'URDD3',LoadCases.LoadFactor;...
                        'URDD5',0 ;'PITCH',0};
                end

            else
                dof = NaN;
                spc = 123456;

                if ~isnan([LoadCases.AileronAngle,LoadCases.ElevatorAngle,LoadCases.RudderAngle])

                    ailer_surf=controlsurfs(cellfun(@(x)contains(x,'aile'),controlsurfs));
                    elev_surf=controlsurfs(cellfun(@(x)contains(x,'elev'),controlsurfs));
                    ruder_surf=controlsurfs(cellfun(@(x)contains(x,'rudder'),controlsurfs));

                    AESTAT = {'ANGLEA',p.Results.AoA;...
%                         'URDD3',LoadCases.LoadFactor;...
%                         'URDD4',0; 'URDD5',0; 'URDD6',0;...
%                         ailer_surf{1}, LoadCases.AileronAngle;...
                        elev_surf{1}, LoadCases.ElevatorAngle};
%                         ruder_surf{1}, LoadCases.RudderAngle;...
%                         'Roll', []; 'Pitch', []; 'Yaw', []};

                else


                    AESTAT = {'ANGLEA',p.Results.AoA;...
                        'URDD3',LoadCases.LoadFactor};

                end

            end
            
            %Return partial data
            SPC_id = ID0;
            TrimData.RefGrid = RefGrid;
            TrimData.SPC_id  = SPC_id;
            TrimData.SPCdof  = spc;
            TrimData.SUPdof  = dof;
            TrimData.Trim    = [];
            
            if isempty(LoadCases)
                return
            end
            
            %Calculate the flight point data - Vectorised
            FlightPoint = getFlightPointData(LoadCases);
            
            %Calculate dynamic pressure
            vel = getAcVelocity(LoadCases);
            q=num2cell(FlightPoint.DynPressure);
            
            %Make the trim entries
            Trim  = repmat(struct('SID', [], 'MACH', [], 'Q', [], 'TrimVar', {{}}), size(LoadCases));
            nTrim = numel(Trim);
            m = {LoadCases.Mach};
            [Trim.MACH] = deal(m{:});
            [Trim.Q]    = deal(q{:});
            for iT = 1 : nTrim
                Trim(iT).TrimVar = [AESTAT',staticControlSurf];
            end
            
            %Set ID numbers for Bulk Data to be called from Case Control
            TRIM_id    = num2cell(ID0 + 1 : SPC_id + nTrim);
            [Trim.SID] = deal(TRIM_id{:});
            AESTAT_id  = ID0 + nTrim + 1 : SPC_id + nTrim + numel(AESTAT(:,1));
            if isempty(MassCases)
                mass_id = [];
            else
                mass_id    = horzcat(MassCases.ID);
            end
            
            %Return the data in a single structure
            TrimData.Trim      = Trim;

            %remove cs from the aestat card
            inda=contains(AESTAT(:,1),'aile');
            inde=contains(AESTAT(:,1),'elev');
            indr=contains(AESTAT(:,1),'rudder');

            ind=inda | inde | indr;
            AESTAT(ind,:)=[];

            TrimData.AESTAT    = AESTAT(:,1)';
            TrimData.AESTAT_id = AESTAT_id;
            TrimData.Mass_id   = mass_id;

            %Check GRAV inout 
            if ~isempty(LoadCases.GRAV)
                TrimData.GRAV=LoadCases.GRAV;
            else
                return
            end
            
        end
    end
    
    methods (Static) %Dynamic Aeroelastic methods
        function FlutterData = getFlutterData(FlightPoints, TrimData, ID0, nModes)
            %getFlutterData Returns a MATLAB structure describing the basic
            %flutter data for the flight points described by the loadcases
            %in 'LoadCases'.
            %
            % The following unsteady aero data is required for a flutter
            % analysis
            %   - Boundary conditions (SPC)
            %   - Eigenanlysis method (EIGRL)
            %   - Unsteady aero data (AERO & MKAERO)
            %   - Flutter data (FLUTTER & FLFACT)
            
            if isa(FlightPoints, 'awi.model.LoadCase')
                FlightPoints = getFlightPointData(FlightPoints);
            else
                validateattributes(FlightPoints, {'awi.model.FlightPoint'}, {'row'}, 'getFlutterData', 'FlightPoints');
            end
            if nargin < 3
                ID0 = 250;
            end
            if nargin < 4
                nModes = [];
            end
            
            %SPC
            SPC.SID   = ID0;
            SPC.C     = TrimData.SPCdof;
            SPC.G     = TrimData.RefGrid(40).ID;
            
            %EIGRL
            EIGRL.SID = SPC.SID + 1;
            EIGRL.V0  = 0;
            EIGRL.V1  = [];
            EIGRL.ND  = 40;
            
            %FLFACT
            %   - Density
            FLFACT(1).SID = EIGRL.SID + 1;
            FLFACT(1).F   = [FlightPoints.DensityRatio];
            FLFACT(1).LAB = 'DENSITY';
            %   - Mach
            FLFACT(2).SID = FLFACT(1).SID + 1;
            FLFACT(2).F   = [FlightPoints.Mach];
            FLFACT(2).LAB = 'MACH';
            %   - Velocity
            FLFACT(3).SID = FLFACT(2).SID + 1;
%             FLFACT(3).F   = [FlightPoints.AcVelocity];

%             FLFACT(3).F   = unique([FlightPoints.AcVelocity, 100 : 20 : 500]);
%             FLFACT(3).F   = 150 : 5 : 450;
            FLFACT(3).F   = 100 : 5 : 400;

%             FLFACT(3).F   = FlightPoints.AcVelocity;
%             FLFACT(3).F   = [FlightPoints.AcVelocity];
            FLFACT(3).LAB = 'VELOCITY';
            
            %FLUTTER
            FLUTTER.SID    = FLFACT(end).SID + 1;
%             FLUTTER.METHOD = 'PKNL';
            % linear analysis 
            FLUTTER.METHOD = 'PK';
            FLUTTER.DENS   = FLFACT(1).SID;
            FLUTTER.MACH   = FLFACT(2).SID;
            FLUTTER.RFREQ  = FLFACT(3).SID;
            
            %TABDMP
            TABDMP1.TID = FLUTTER.SID(end) + 1;
            TABDMP1.x   = [0   , 1000];
            TABDMP1.y   = [0.02, 0.02];
            
            %MKAERO
            MKAERO.M = unique([FlightPoints.Mach]);
            MKAERO.K = [0.001, 0.002, 0.005, 0.01, 0.02, 0.03, 0.06, 0.1, 0.2, 0.4, 0.7, 1.0, 1.3, 1.6, 2, 2.4, 3];
%             MKAERO.K = [0.001, 0.02 , 0.04, 0.06, 0.08, 0.1, 0.12, 0.2, 0.5, 1, 1.5, 2, 2.5, 3, 3.5];
            %MKAERO.K = [0.005, 0.001, 0.01 : 0.01 : 3.5]; 
            
            %SET1 (For tracking flutter modes in a SOL 200 analysis)
            if isempty(nModes)
                SET1 = [];
            else
                nRigidBidyDOF = 6 - numel(SPC.C);
                SET1.SID = TABDMP1.TID + 1;
                SET1.ID  = [nRigidBidyDOF + 1, nRigidBidyDOF + nModes];
            end
            
            %Collect
            FlutterData.SPC     = SPC;
            FlutterData.EIGRL   = EIGRL;
            FlutterData.FLFACT  = FLFACT;
            FlutterData.FLUTTER = FLUTTER;
            FlutterData.MKAERO  = MKAERO;
            FlutterData.TABDMP1 = TABDMP1;
            FlutterData.SET1    = SET1;
            FlutterData.Mass_id = [];
            
        end
        function GustData = getGustData(TrimData,FlightPoints,LoadCases)
            %getGustData Returns a MATLAB structure describing the basic
            %gust data for the flight points described by the loadcases
            %in 'LoadCases'.
            %
            % The following unsteady aero data is required for a gust
            % analysis
            %   - Boundary conditions (SPC)
            %   - Eigenanlysis method (EIGRL)
            %   - Unsteady aero data (AERO & MKAERO)
            

             % determine 1MC gust
             L=LoadCases.GustLength; % gust length
             N=numel(L);
             Ng=64; %number of points
             bulk_data=zeros(N,Ng*2);
             Uds=zeros(1,N); % data record for gust amplitude
             
            % loop for gust families: positive only -------------------
             for n=1:N
                 H=L(n)*0.5;
                 Kg=1; % gust alleviation factor
                 
                 tmax=L(n)/FlightPoints.AcVelocity;
                 ts=linspace(0,tmax,Ng);
                 %Calculate gust velocity
                 Uref = [ ...
                     0    , 15000, 60000 ; ...
                     17.07, 13.41, 6.36 ];
                 Ug = interp1(Uref(1, :), Uref(2, :),FlightPoints.Altitude);  %[m/s], EAS
                 U_ref = Ug ./ sqrt(FlightPoints.DensityRatio);  %[m/s], TAS
                 
                 Uds(n)=U_ref*Kg*(H/107)^(1/6);
                 
                 gust_1mc=Uds(n)*(1-cos(2*pi*FlightPoints.AcVelocity*ts/L(n)))*0.5;                 
                 
                 for i=0:numel(ts)-1
                     
                     index1=1+2*i;
                     index2=2+2*i;
                     bulk_data(n,index1)=ts(i+1);
                     bulk_data(n,index2)=gust_1mc(i+1)/max(abs(gust_1mc));
                     
                 end
                 
             end           
             GustData.bulk_data=bulk_data;
             GustData.amplitudes=LoadCases.GustDirection*Uds;
             
             % enter TLOAD data
             GustData.TLid=100;
             
             % TABLED1 ID
             GustData.TABLED1.id=1;
             
             % enter DLOAD data
             GustData.DLid=10;
             GustData.DLtab=1;
             GustData.DLS=1.0; % scale factor 1
             GustData.DLSi=1.0; % scale factor 2
             
             % enter DAREA data
             GustData.DAREA.sid=200;
             % TODO - this is not very generic (was 34)...
             GustData.DAREA.p1=TrimData.RefGrid(1).GID; % choose the wing root Grid ID
             GustData.DAREA.c1=1;
             GustData.DAREA.a1=1;
             
             % enter GUST data
             GustData.GUSTid=300;
             GustData.Wg=LoadCases.GustDirection.*Uds/FlightPoints.AcVelocity;
             GustData.X0=-100; %initial position of the ac
             
             % FREQUQ
             GustData.FREQ.id=600;
             GustData.FREQ.F1=0;
             GustData.FREQ.DF=0.05;
             GustData.FREQ.NDF=600;
             
             % TSTEP
             GustData.TSTEP.id=700;
             GustData.TSTEP.N=200;
             GustData.TSTEP.DT=0.025;
             GustData.TSTEP.NO1=1;
             
             % EIGR
             GustData.EIGRL.SID = 20;
             GustData.EIGRL.V0  = 0;
             GustData.EIGRL.V1  = [];
             GustData.EIGRL.ND  = 30;
             
             % TABDAMP1           
             GustData.TABDMP1.TID = 1000;
             GustData.TABDMP1.x   = [0   , 1000];
             GustData.TABDMP1.y   = [0.12, 0.12];
             
             dat = [GustData.TABDMP1.x ;GustData.TABDMP1.y];
             str = cellstr(num2str(dat(:), '%#-8.2g'));
             GustData.str = [str ; {'ENDT'}];
             if numel(str) > 8
                 error('Update code for writing structural damping table');
             end
        end
        function writeFlutterEntries(fid, FlutterData)
            %writeFlutterEntries
            awi.fe.FEBaseClass.writeColumnDelimiter(fid, 'normal');
            %spc
            comp = FlutterData.SPC.C;
            if ischar(comp)
                comp = str2double(comp);
            end
            % split line
            awi.fe.FEBaseClass.writeColumnDelimiter(fid, 'normal');
            
            % method - EIGR (charles)
            if isempty(FlutterData.EIGRL.V1)
                fprintf(fid, '%-8s%-8i%-8s%-#8.3g%-16s%-8i\r\n', 'EIGR', FlutterData.EIGRL.SID, 'MGIV', 0, blanks(16),FlutterData.EIGRL.ND);
            else
                fprintf(fid, '%-8s%-8i%-8s%-#8.3g%-#8.3g\r\n', 'EIGR', FlutterData.EIGRL.SID, 'MGIV',0, FlutterData.EIGRL.V1); %, FlutterData.EIGRL.ND);
            end
            %flutter
            fprintf(fid, '%-8s%-8i%-8s%-8i%-8i%-8i%-8s%-8s%-8s\r\n', 'FLUTTER', ...
                FlutterData.FLUTTER.SID  , FlutterData.FLUTTER.METHOD, ...
                FlutterData.FLUTTER.DENS , FlutterData.FLUTTER.MACH  , ...
                FlutterData.FLUTTER.RFREQ, 'L', blanks(8), '1.E-4');
            %flfact
            awi.methods.Nastran.writeFLFACT(fid, FlutterData.FLFACT);
            %mkaero
            awi.methods.Nastran.writeMKAERO1(fid, FlutterData.MKAERO);
            %set1
            if isfield(FlutterData, 'SET1') && ~isempty(FlutterData.SET1)
                if isfield(FlutterData.SET1, 'SpecificModes')
                    fprintf(fid, '%-8s%-8i%-8i%-8i\r\n', 'SET1', FlutterData.SET1.SID, ...
                        FlutterData.SET1.ID(1), FlutterData.SET1.ID(2));
                else
                    fprintf(fid, '%-8s%-8i%-8i%-8s%-8i\r\n', 'SET1', FlutterData.SET1.SID, ...
                        FlutterData.SET1.ID(1), 'THRU', FlutterData.SET1.ID(end));
                end                
            end
            %tabdmp1
            dat = [FlutterData.TABDMP1.x ; FlutterData.TABDMP1.y];
            str = cellstr(num2str(dat(:), '%#-8.2g'));
            str = [str ; {'ENDT'}];
            if numel(str) > 8
                error('Update code for writing structural damping table');
            end
            fprintf(fid, ['%-8s%-8i%-8s\r\n%-8s', repmat('%-8s', [1, numel(str)]), '\r\n'], ...
                'TABDMP1', FlutterData.TABDMP1.TID, 'G', blanks(8), str{:});
            
        end
        function writeUnsteadyAeroEntry(fid, Aircraft, LoadCases, varargin)
            %writeUnsteadyAeroEntry Writes the equivalent MSC.Nastran AERO
            %entry using the data in 'LoadCase'.
            %
            %   - Assumes that the aero coordinate system is the basic
            %     coordinate system.
            
            p = inputParser;
            addParameter(p, 'SymXZ', false, @(x)validateattributes(x, {'logical'}, {'scalar'}));
            addParameter(p, 'SymXY', false, @(x)validateattributes(x, {'logical'}, {'scalar'}));
            parse(p, varargin{:});
            
            %Write the data to the file
            fprintf(fid, '%-8s%-8i%-8.4f%-8.4f%-8.4f%-8i%-8i\r\n', ...
                'AERO', 0, LoadCases(1).AcVelocity, Aircraft.RefChord, ...
                LoadCases(1).RefDensity, p.Results.SymXZ, p.Results.SymXY);
            
        end
    end
    
    methods (Static) %Generic optimisation methods
        function DesignModel = extractDesignModel(pchFile)
            %extractDesignModel Extracts the bulk data from the PUNCH file
            %(.pch) which is output by SOL200 during the optimisation run.
            %
            % Syntax :
            %   DesignModel = extractDesignModel(pchFile)
            %
            % Inputs :
            %   - Required :
            %       + 'pchFile' - Name of the PUNCH file containing the
            %                     design model data.
            %
            % Outputs :
            %   - 'DesignModel' - Structure containing the bulk data
            %                     information.
            %
            % Author    : Christopher Szczyglowski
            % Username  : cs1414
            % Email     : cs1414@bristol.ac.uk
            % Timestamp : 20-May-2019 11:55:04
            %
            % Copyright (c) 2019 Christopher Szczyglowski            % All Rights Reserved
            %
            %
            % Revision: 1.0 20-May-2019 11:55:04
            %   - Initial function :
            %
            
            [~, ~, ext] = fileparts(pchFile);
            validatestring(ext, {'.pch', '.txt', '.dat', '.bdf'}, 'extractDesignModel', 'file extension');
            assert(exist(pchFile, 'file') == 2, sprintf(['File ''%s'' ', ...
                'not found. Unable to continue.'], pchFile));
            
            %% Import data from the .pch file
            
            %Read every line from the file (assume all in one block) & discard comments
            fid = fopen(pchFile, 'r');
            C = textscan(fid, '%s', 'delimiter', '\n','CommentStyle', '$', 'whitespace', '');
            bulk = C{1};
            clear C
            fclose(fid);
            
            %Determine unique bulk data names
            col1           = cellfun(@(x) x(1 : 8), bulk, 'Unif', false);
            bulkNames      = col1(~startsWith(col1, '*'));
            idx            = contains(bulkNames, '*');
            bulkNames(idx) = cellfun(@(x) x(1 : strfind(x, '*') - 1), bulkNames(idx), 'Unif', false);
            bulkNames      = strtrim(unique(bulkNames));
            
            %Grab formatting data
            FormatData = i_grabBulkFormat;
            
            %% Convert raw text data into MSC.Nastran bulk data
            %   - Assume large-field format
            
            for iB = 1 : numel(bulkNames)
                
                %Grab formatting data for this card
                bName  = bulkNames{iB};
                if ~isfield(FormatData, bName)
                    warning(['New bulk data type %s has been skipped. ', ...
                        'Add to the ''BulkFormat'' variable.'], bName);
                    continue
                end
                
                %Find all lines containing data
                idxBulk = startsWith(bulk, bName);
                ind1  = find(idxBulk, 1, 'first');
                indN  = find(idxBulk, 1, 'last');
                bulkIndex = ind1 : indN + numel(FormatData.(bName).Fields) - 1;
                bd    = bulk(bulkIndex);
                
                %Extract the data
                DesignModel.(bName) = i_extractBulkData(bd, FormatData.(bName));
                
            end
            
            function FormatData = i_grabBulkFormat
                
                %% DESVAR
                DESVAR.Fields = {{'ID', 'LABEL', 'XINIT', 'XLB'} ; {'XUB', 'DELXV'}};
                DESVAR.Format = {'rcrr' ; 'rr'};
                FormatData.DESVAR = DESVAR;
                
                %% PBEAM
                %   - Need to check if "SO" is YES or YESA as this affects how many entries
                %   are printed in the punch file.
                % PBEAM.Fields = { ...
                %     {'PID' , 'MID' , 'A_A' , 'I1_A'} ; {'I2_A' , 'I12_A', 'J_A' , 'NSM_A'} ; ...
                %     {'C1_A', 'C2_A', 'D1_A', 'D2_A'} ; {'E1_A' , 'E2_A' , 'F1_A', 'F2_A'}  ; ...
                %     {'SO'  , 'X_XB', 'A_B' , 'I1_B'} ; {'I2_B' , 'I12_B', 'J_B' , 'NSM_B'} ; ...
                %     {'C1_B', 'C2_B', 'D1_B', 'D2_B'} ; {'E1_B' , 'E2_B' , 'F1_B', 'F2_B'}  ; ...
                %     {'K1'  , 'K2'  , 'S1'  ,'S2'   } ; {'NSI_A', 'NSI_B', 'CW_A', 'CW_B'}};
                PBEAM.Fields = { ...
                    {'PID' , 'MID' , 'A_A' , 'I1_A'} ; {'I2_A' , 'I12_A', 'J_A' , 'NSM_A'} ; ...
                    {'C1_A', 'C2_A', 'D1_A', 'D2_A'} ; {'E1_A' , 'E2_A' , 'F1_A', 'F2_A'}  ; ...
                    {'SO'  , 'X_XB', 'A_B' , 'I1_B'} ; {'I2_B' , 'I12_B', 'J_B' , 'NSM_B'} ; ...
                    {'K1'  , 'K2'  , 'S1'  ,'S2'   } ; {'NSI_A', 'NSI_B', 'CW_A', 'CW_B'}  ; ...
                    {'M1_A', 'M2_A', 'M1_B', 'M2_B'} ; {'N1_A' , 'N2_A' , 'N1_B', 'N2_B'} };
                PBEAM.Format = { ...
                    'rrrr' ; 'rrrr' ; 'rrrr' ; 'rrrr' ; 'crrr' ; 'rrrr' ; 'rrrr' ; 'rrrr' ; 'rrrr' ; 'rrrr' };
                FormatData.PBEAM = PBEAM;
                
                %% PDAMP
                PDAMP.Fields = {{'PID', 'B'}};
                PDAMP.Format = {'rr'};
                FormatData.PDAMP = PDAMP;
                
            end
            
            function BulkData = i_extractBulkData(rawBulkText, FormatData)
                %extractBulkData Extracts the data from the raw text and
                %returns a MATLAB structure.
                
                fNames = FormatData.Fields;
                format = FormatData.Format;
                
                %Extract data from each line (vectorised)
                nData  = numel(format);
                nLines = numel(rawBulkText);
                
                %How many full lines?
                nCol  = cellfun(@numel, fNames);
                idxF  = nCol == 4;
                nFull = nnz(idxF);
                nPart = find(~idxF);
                if numel(nPart) > 1
                    error('Update code to loop through all partial lines');
                end
                
                %% Extract data for the full-lines
                %Index the full lines only
                fIndex = arrayfun(@(i) i : nData : nLines, 1 : nFull, 'Unif', false);
                fIndex = vertcat(fIndex{:});
                fIndex = sort(fIndex(:));
                fData  = rawBulkText(fIndex);
                
                %Trim columns 1 & 10
                fData = cellfun(@(x) x(9 : 72), fData, 'Unif', false);
                
                %Split each line into 4 columns
                fData = cellfun(@(x) {x(1 : 16), x(17 : 32), x(33 : 48) x(49 : 64)}, fData, 'Unif', false);
                fData = vertcat(fData{:});
                
                %Convert data to correct type
                fNam = horzcat(fNames{idxF});
                form = horzcat(format{idxF});
                data = formatRawData(fNam, form, fData);
                
                %% Extract data for the partial lines
                if isempty(nPart)
                    data_ = [];
                else
                    pIndex = nPart : nData : nLines;
                    pData  = rawBulkText(pIndex);
                    
                    %Trim column 1
                    pData = cellfun(@(x) x(9 : end), pData, 'Unif', false);
                    
                    %Define bounds for indexing - Assume full data provided
                    nCol  = numel(format{~idxF});
                    ub    = cumsum(repmat(16, [1, nCol]));
                    lb    = [1, ub(1 : end - 1) + 1];
                    
                    %How much data is actually provided? Pad with blank data so we can
                    %still vectorise
                    nCol   = cellfun(@numel, pData);
                    maxCol = max(nCol);
                    del    = (maxCol - nCol) ./ 16;
                    udel   = unique(del);
                    udel   = udel(udel ~= 0);
                    for i = 1 : numel(udel)
                        idx = del == udel(i);
                        pData(idx) = strcat(pData(idx), {repmat(blanks(16), [1, udel(i)])});
                    end
                    
                    %Check that we have extracted the right amount of data
                    lb = lb(lb < maxCol);
                    ub = ub(1 : numel(lb));
                    
                    %Split each line into columns
                    pData = arrayfun(@(i) cellfun(@(x) x(lb(i) : ub(i)), pData, 'Unif', false), 1 : numel(lb), 'Unif', false);
                    pData = horzcat(pData{:});
                    
                    %Convert data to correct type
                    pNam  = horzcat(fNames{~idxF});
                    form  = horzcat(format{~idxF});
                    form  = form(1 : numel(lb));
                    data_ = formatRawData(pNam, form, pData);
                end
                
                %% Collect data from full & partial lines
                data     = [data, data_];
                bNames   = horzcat(fNames{:});
                BulkData = cell2struct(data, bNames, 2);
                
                function data = formatRawData(name, format, rawdata)
                    %formatRawData Returns the cell data 'rawdata' as a collection of
                    %either numeric or cell-string data.
                    
                    if isempty(rawdata)
                        data = [];
                        return
                    end
                    
                    %Modify size of 'rawdata' to group together data from each card
                    nSet = numel(name) / 4;
                    if nSet > 1
                        ub_   = nSet : nSet : size(rawdata, 1);
                        lb_   = [1, ub_(1 : end - 1) + 1];
                        rawdata = arrayfun(@(i) rawdata(lb_(i) : ub_(i), :)', 1 : numel(lb_), 'Unif', false);
                        rawdata = cellfun(@(x) x(:)', rawdata, 'Unif', false);
                        rawdata = vertcat(rawdata{:});
                    end
                    
                    %Look for empty data
                    idx = cellfun(@(x) isequal(x, blanks(16)), rawdata);
                    rawdata(idx) = {pad('NaN', 16)};
                    
                    %Format real & string data
                    data = cell(size(name));
                    idxR = format == 'r';
                    data(idxR)  = arrayfun(@(i) str2num(vertcat(rawdata{:, i}))', find(idxR), 'Unif', false); %#ok<*ST2NM>
                    data(~idxR) = arrayfun(@(i) rawdata(:, i)', find(~idxR), 'Unif', false);
                    
                end
                
            end
            
        end
    end
    
    methods (Static) %Beam Optimisation methods
        function printTrimSubcases(fid, sc0, TrimData, DesConstr, pad, fileNo)
            %printTrimSubcases Prints each trim subcase into the file with
            %identifier 'fid'.
            
            scInd = sc0 : sc0 + numel(TrimData.Trim) - 1;
            
            for ii =  1 : numel(scInd)
                fprintf(fid, '$S U B C A S E  %i\r\n', scInd(ii));
                fprintf(fid, 'SUBCASE %i\r\n'        , scInd(ii));
                fprintf(fid, '%sANALYSIS       = SAERO\r\n', pad);
                fprintf(fid, '%sTRIM           = %i\r\n'   , pad, TrimData.Trim(ii).SID);
                fprintf(fid, '%sSPC            = %i\r\n'   , pad, TrimData.SPC_id);
                fprintf(fid, '%sDESSUB         = %i\r\n'   , pad, DesConstr.ID);
                fprintf(fid, '%sFORCE(NOPRINT) = ALL\r\n'  , pad);
                fprintf(fid, '%sDISP(NOPRINT)  = ALL\r\n'  , pad);
                if ~isempty(fileNo)
                    fprintf(fid, '%sTRIMF(UNIT=%i) = ALL\r\n'  , pad, fileNo(ii));
                end
                if ~isempty(TrimData.Mass_id)
                    fprintf(fid, '%sMASSSET = %i\r\n', pad, TrimData.Mass_id(ii));
                end
            end
            
        end
        function printStaticSubcases(fid, sc0, StaticData, DesConstr, pad)
            %printTrimSubcases Prints each static subcase into the file
            %with identifier 'fid'.
            
            scInd = sc0 : sc0 + numel(StaticData.LOAD_id) - 1;
            
            for ii = 1 : numel(scInd)
                fprintf(fid, '$S U B C A S E  %i\r\n', scInd(ii));
                fprintf(fid, 'SUBCASE %i\r\n'        , scInd(ii));
                fprintf(fid, '%sANALYSIS       = STATICS\r\n', pad);
                fprintf(fid, '%sLOAD           = %i\r\n'     , pad, StaticData.LOAD_id(ii));
                fprintf(fid, '%sSPC            = %i\r\n'     , pad, StaticData.SPC_id(ii));
                fprintf(fid, '%sDESSUB         = %i\r\n'     , pad, DesConstr.ID);
                fprintf(fid, '%sFORCE(NOPRINT) = ALL\r\n'    , pad);
                fprintf(fid, '%sDISP(NOPRINT)  = ALL\r\n'    , pad);
                if ~isempty(StaticData.Mass_id)
                    fprintf(fid, '%sMASSSET = %i\r\n', pad, StaticData.Mass_id(ii));
                end
                
            end
        end
        function printBucklingSubcases(fid, TrimData, pad, BuckData, GBConstr)
            %printTrimSubcases Prints each buckling subcase into the file
            %with identifier 'fid'.
            %
            %   - Each buckling analysis requires two subcases
            %       1. Static analysis using the applied trim loads from a
            %       previous subcase to determine the differential
            %       stiffness matrix.
            %       2. Linear buckling analysis to calculate the buckling
            %          eigenvalue and the necessary constraints.
            
            nTr  = numel(TrimData.Trim);
            scNo = nTr+1 : nTr + nTr*2;
            scNo = reshape(scNo, [2, numel(scNo)/2]);
            for ii = 1 : numel(TrimData.Trim)
                %SOL 101 - Apply the trim loads to get differential stiffness matrix
                fprintf(fid, '$S U B C A S E  %i\r\n', scNo(1, ii));
                fprintf(fid, 'SUBCASE %i\r\n'        , scNo(1, ii));
                fprintf(fid, '%sANALYSIS          = STATICS\r\n', pad);
                fprintf(fid, '%sLOAD              = %i\r\n'     , pad, ii);
                fprintf(fid, '%sSPC               = %i\r\n'     , pad, BuckData.SPC_id);
                if ~isempty(BuckData.Mass_id)
                    fprintf(fid, '%sMASSSET = %i\r\n', pad, BuckData.Mass_id(ii));
                end
                %SOL 105 - Buckling accounting for differential stiffness
                fprintf(fid, '$S U B C A S E  %i\r\n', scNo(2, ii));
                fprintf(fid, 'SUBCASE %i\r\n'        , scNo(2, ii));
                fprintf(fid, '%sANALYSIS          = BUCK\r\n', pad);
                fprintf(fid, '%sSPC               = %i\r\n'  , pad, BuckData.SPC_id);
                fprintf(fid, '%sMETHOD            = %i\r\n'  , pad, BuckData.EIGRL_id);
                fprintf(fid, '%sSTATSUB(BUCKLING) = %i\r\n'  , pad, scNo(1, ii));
                fprintf(fid, '%sDISP(NOPRINT)     = ALL\r\n' , pad);
                fprintf(fid, '%sDESSUB            = %i\r\n'  , pad, GBConstr.ID);
                if ~isempty(BuckData.Mass_id)
                    fprintf(fid, '%sMASSSET = %i\r\n', pad, BuckData.Mass_id(ii));
                end
            end
            
        end
        function printFlutterSubcases(fid, FlutterData, FConstr, scNo, pad, bDisp)
            %printFlutterSubcases Prints each flutter subcase into the file
            %with identifier 'fid'.
            
            if nargin < 6
                bDisp = false;
            end
            
            fprintf(fid, '$S U B C A S E  %i\r\n', scNo);
            fprintf(fid, 'SUBCASE %i\r\n'        , scNo);
            if ~isempty(FConstr) %SOL200 case control 
                fprintf(fid, '%sANALYSIS = FLUTTER\r\n', pad);
                fprintf(fid, '%sDESSUB   = %i\r\n'  , pad, FConstr.ID);
            end
            
            if contains(FlutterData.SPC.C,'123456')
                fprintf(fid, '%sSPC      = %i\r\n'  , pad, FlutterData.SPC.SID);
            end
            
            
            fprintf(fid, '%sMETHOD   = %i\r\n'  , pad, FlutterData.EIGRL.SID);
            fprintf(fid, '%sFMETHOD  = %i\r\n'  , pad, FlutterData.FLUTTER.SID);
            fprintf(fid, '%sSDAMP    = %i\r\n'  , pad, FlutterData.TABDMP1.TID);
            if bDisp %Request flutter modeshapes
                fprintf(fid, '%sDISP     = ALL\r\n' , pad);
            end
            if ~isempty(FlutterData.Mass_id) %Mass sets
                fprintf(fid, '%sMASSSET = %i\r\n', pad, FlutterData.Mass_id);
            end
        end
        function printDesOptParam(fid, Options)
            %printDesOptParam Prints the data in 'DesOptParam' to the file
            %with identifier 'fid'.
            
            %How many parameters
            nP = size(Options.Sol200ParamMap, 2);
            
            DesOptParam = [ ...
                Options.Sol200ParamMap(2, :) ; get(Options, Options.Sol200ParamMap(1, :))]';
            
            %Convert numeric data to strings...
            %   - Deal with integers
            idx = ismember(DesOptParam(:, 1), {'DESMAX', 'NASPR0'});
            DesOptParam(idx, 2) = cellstr(num2str(vertcat(DesOptParam{idx, 2}), '%-8i'));
            %   - Parameter 'DELOBJ' doesn't seem to work well when
            %     printed using %g, use %f instead
            idx_do = ismember(DesOptParam(:, 1), 'DELOBJ');
            DesOptParam(idx_do, 2) = cellstr(num2str(vertcat(DesOptParam{idx_do, 2}), '%-8f'));
            %   - Everything else is simple...
            idx = not(or(idx, idx_do));
            DesOptParam(idx, 2)  = cellstr(num2str(vertcat(DesOptParam{idx, 2}), '%#-8.3g'));
            
            %Get rid of DABOBJ --> Messes up optimisation run
            idx = ismember(DesOptParam(:, 1), 'DABOBJ');
            DesOptParam = DesOptParam(~idx, :);
            
            %Determine number of full lines
            if mod(nP, 4) ~= 0
                nLines = floor(nP ./ 4);
                rem    = DesOptParam(nLines * 4 + 1 : end, :)';
                lineData = DesOptParam(1 : nLines * 4, :)';
            else
                nLines   = floor(nP ./ 4);
                if nLines == 0
                    nLines = 1;
                end
                lineData = DesOptParam';
                rem = {};
            end
            
            %Reshape to prepare for vectorised printing
            lineData = lineData(:);
            lineData = reshape(lineData, [8, nLines]);
            
            %Add column 1 data
            lineData = [{'DOPTPRM'}, repmat({blanks(8)}, [1, nLines - 1]) ; lineData];
            
            %Define format
            fmt = ['%-8s', repmat('%-8s', [1, 8]), '\r\n'];
            
            %Print
            fprintf(fid, fmt, lineData{:});
            if ~isempty(rem)
                fmt = ['%-8s', repmat('%-8s%-8s', [1, size(rem, 2)]), '\r\n'];
                fprintf(fid, fmt, blanks(8), rem{:});
            end
            
        end
        function [TrimData, StaticData, BuckData, FlutterData] = getAnalysisData(FEM, Aircraft, LoadCases, MassCases, Options)
            %getAnalysisData Returns a set of MATLAB structures containing
            %the data required to write the Bulk Data and Case Control
            %entries for the trim, buckling and flutter cases.
            
            %Initial ID of extra bulk data entries
            ID0 = 250;
            
            %Where will constraints be applied?
            AllFEM = flatlist(FEM);
            Fuselage = AllFEM(ismember({AllFEM.Name}, 'Fuselage'));
            if ~isempty(Fuselage)
                RefGrid = Fuselage.Nodes(1);
            else
                error(['Unable to proceeed as there is no reference ', ...
                    'node for the constraints to be applied to.']);
            end
            
            type  = get(LoadCases, {'LoadCaseType'});
            idxS  = ismember(type, 'Static');
            idxPG = ismember(type, 'Pratt Gust');
            idxM  = ismember(type, {'Manoeuvre', 'Pratt Gust'});
            
            %Update load factor for Pratt Gust
            calculateGustLoadFactor(LoadCases(idxPG), Aircraft);
            
            %Get the TrimData for the static aeroelastic load cases
            TrimLoadCases = LoadCases(idxM);
            if isempty(MassCases)
                TrimMassCases = [];
            else
                TrimMassCases = MassCases(idxM);
            end
            TrimData = awi.methods.Nastran.getTrimData(FEM, Aircraft, TrimLoadCases, TrimMassCases, ID0, RefGrid);
            
            ID0 = ID0 + nnz(idxM) + 1;
            
            %Define the StaticData for the static load cases
            if any(idxS)
                StaticData.RefGrid = TrimData.RefGrid;
                StaticData.SPC_id  = ID0 + 1;
                StaticData.SPCdof  = TrimData.SPCdof;
                StaticData.GRAV_id = StaticData.SPC_id + 1;
                StaticData.LOAD_id = StaticData.GRAV_id + 1 : StaticData.GRAV_id + nnz(idxS);
                if isempty(MassCases)
                    StaticData.Mass_id = [];
                else
                    StaticData.Mass_id = horzcat(MassCases(idxS).ID);
                end
                ID0 = StaticData.LOAD_id(end) + 1;
            else
                StaticData = [];
            end
            
            %Define the BuckData for the linear buckling load cases
            if Options.GlobalBuckConstr
                BuckData.RefGrid  = TrimData.RefGrid;
                BuckData.SPC_id   = ID0 + 1;
                BuckData.SPCdof   = TrimData.SPCdof;
                BuckData.EIGRL_id = BuckData.SPC_id + 1;
                BuckData.NModes   = 3;
                BuckData.Mass_id  = TrimData.Mass_id;
                ID0 = BuckData.EIGRL_id;
            else
                BuckData = [];
            end
            
            %Get the flutter data for the flutter load cases
            FlutterData = awi.methods.Nastran.getFlutterData(LoadCases(1), TrimData, ID0);
            %   - Structural damping
            FlutterData.TABDMP1.TID = FlutterData.FLUTTER.SID + 1;
            FlutterData.TABDMP1.x   = [0, 50];
            FlutterData.TABDMP1.y   = [0.03, 0.03];
            %   - Define modes
            FlutterData.SET1.SID = FlutterData.TABDMP1.TID + 1;
            if Options.IgnoreRigidBodyModes
                %Calculate how many rigid body DOF there are
                m0 = 6 - numel(num2str(TrimData.SPCdof));
            else
                m0 = 0;
            end
            FlutterData.SET1.ID  = m0 + 1 : m0 + Options.NumFlutterModes;
            if isfield(TrimData, 'Mass_id')
                FlutterData.Mass_id   = TrimData.Mass_id;
            else
                FlutterData.Mass_id = [];
            end
            
        end
    end
    
    methods(Static)
        function controlsurfs = getControlSurfs(FEM)
           %What aero controllers have been defined?
            FEM = flatlist(FEM);
            AeroControls = [FEM.AeroControlSurf];
            if isempty(AeroControls)
                controlsurfs = {};
            else                                
                controlsurfs = {AeroControls.LABEL};                
            end 
        end
    end
    
    
    methods (Static) %GFEM Optimisation methods
        
    end
    
    methods (Static) %Writing MSC.Nastran bulk data
        function writeFLFACT(fid, FLFACT, bComment)
            %writeFLFACT : Writes the entries for a FLFACT card
            
            assert(and(fid ~= -1, numel(fid) == 1), ['''fid'' must be ', ...
                'a valid file identifier.']);
            assert(isstruct(FLFACT), 'Expected FLFACT to be a structure.');
            assert(all(arrayfun(@(s) isfield(s, 'SID'), FLFACT)), ...
                'Each FLFLACT entry must have an ''SID'' field.');
            assert(all(arrayfun(@(s) isfield(s, 'F'), FLFACT)), ...
                'Each FLFLACT entry must have an ''F'' field.');
            assert(all(arrayfun(@(s) isfield(s, 'LAB'), FLFACT)), ...
                'Each FLFLACT entry must have an ''LAB'' field.');
            
            if nargin < 3
                bComment = false;
            end
%             write_text_header( fid, 'FLFACT : Specifies density ratios, Mach numbers, reduced frequencies and velocities for flutter analysis', 'sub')
%             write_column_headings(fid)
            if bComment
                fprintf(fid, '$       SID     F1      F2      F3      F4      F5      F6      F7      F8\r\n');
            end
            
            strFormat = '%-8.3f';

            for i = 1 : numel(FLFACT)
                nEntry    = numel(FLFACT(i).F);             % total no. entries
                nLines    = ceil((nEntry - 7) ./ 8) + 1;    % total no. lines
                if nLines == 1
                    nEntryEnd = nEntry;
                else
                    nEntryEnd = mod((nEntry - 7), 8);           % no. entries on final line
                end
                maxEntry  = [7, repmat(8, [1, nLines-1])];  % max no. entries on each line
                %Define number of entries on each line
                if nLines == 1
                    nDat = nEntry;
                else
                    nDat =[7, repmat(8, [1, nLines-2]), nEntryEnd];
                end
                %Define bounds
                ub = cumsum(nDat);
                lb = [1, ub(1 : end-1) + 1];
                %Loop through lines and write data to file
                for iL = 1 : nLines
                    if iL == nLines
                        %Construct blanks + label string
                        nBlankCol = maxEntry(iL) - nDat(iL);
                        str = sprintf('%s%s',repmat(blanks(8), [1, nBlankCol]), FLFACT(i).LAB);
                    else
                        str = blanks(8);
                    end
                    %Print FLFACT
                    if iL == 1
                        fprintf(fid, '%-8s%-8i%-s%-s\r\n', 'FLFACT', FLFACT(i).SID, ...
                            sprintf(strFormat, FLFACT(i).F(lb(iL) : ub(iL))), str);
                    else
                        fprintf(fid, '%-8s%-s%-s\r\n', blanks(8), ...
                            sprintf(strFormat, FLFACT(i).F(lb(iL) : ub(iL))), str);
                    end
                end
            end
            
        end   
        function writeMKAERO1(fid, MKAERO1, bComment)
            %writeMKAERO1 : Writes the entries for a MKAERO1 card
            
            assert(and(fid ~= -1, numel(fid) == 1), ['''fid'' must be ', ...
                'a valid file identifier.']);
            assert(isstruct(MKAERO1), 'Expected MKAERO1 to be a structure.');
            assert(all(arrayfun(@(s) isfield(s, 'M'), MKAERO1)), ...
                'Each MKAERO1 entry must have an ''M'' field.');
            assert(all(arrayfun(@(s) isfield(s, 'K'), MKAERO1)), ...
                'Each MKAERO1 entry must have an ''K'' field.');
            if nargin < 3
                bComment = false;
            end
            
            nEntry = 8; % number of M/K vlaues per line
            
            k = MKAERO1.K;
            m = MKAERO1.M;
            
            %Determine how many MKAERO1 entries need to be printed            
            %   - Total number of Mach/k values
            nM = numel(m);
            nK = numel(k);            
            %   - Number of full lines required
            nM_full = floor(nM / nEntry);
            nK_full = floor(nK / nEntry);            
            %   - Amount of data on the final (non-full) line
            rem_M = mod(nM, nEntry);
            rem_K = mod(nK, nEntry);            
            %   - If the remainder is zero then change to 1 and accept that
            %     there will be duplicate data
            if rem_M == 0
                rem_M = 1;
                m(end+1) = m(1)*1.0;
            end
            if rem_K == 0
                rem_K = 1;
                k(end+1) = k(end)*1.0;
            end            
            %	- Total number of cards
            nCard = ceil(max([nM, nK]) ./ nEntry);
            
            %Define upper and lower bounds for indexing the data
            nCol_M   = [repmat(nEntry, [nM_full, 1]) ; rem_M];
            nCol_K   = [repmat(nEntry, [nK_full, 1]) ; rem_K];
            bounds_M = [cumsum(nCol_M)-nCol_M + 1 , cumsum(nCol_M)];
            bounds_K = [cumsum(nCol_K)-nCol_K + 1 , cumsum(nCol_K)];
            
            %Print data to file
            
            % MKAERO1 --> Table defining mach numbers and reduced frequency values for
            % calculation of aerodynamic matricies
            if bComment
                comment = 'MKAERO1 : Mach Number - Frequency Table';
%                 writeTextHeader(fid, comment, 'sub')
%             writeColumnHeadings(fid);
%             %
%             fprintf(fid, '$       m1      m2      m3      m4      m5      m6      m7      m8\r\n');
%             fprintf(fid, '$       k1      k2      k3      k4      k5      k6      k7      k8\r\n');
            end            
            
            for iC = 1 : nCard
                if iC <= nM_full + 1
                    fprintf(fid, '%-8s%s\r\n', 'MKAERO1', num2str(m(bounds_M(iC, 1) : bounds_M(iC, 2)), '%-8.5f'));
                else
                    % pad with the first value of M
                    %         fprintf(fid, '%-8s%s\r\n', 'MKAERO1', num2str(MKAERO1.M(end), '%-8.5f'));
                    fprintf(fid, '%-8s%s\r\n', 'MKAERO1', num2str(m(end)*1.0, '%-8.5f')); %avoid duplicates
                end
                if iC <= nK_full + 1
                    fprintf(fid, '%-8s%s\r\n', blanks(8), num2str(k(bounds_K(iC, 1) : bounds_K(iC, 2)), '%#-8.5f'));
                else
                    % pad with the first value of k
                    %         fprintf(fid, '%-8s%s\r\n', 'MKAERO1', num2str(MKAERO1.k(end), '%#-8.5f'));
                    fprintf(fid, '%-8s%s\r\n', 'MKAERO1', num2str(k(end)*1.01, '%#-8.5f'));
                end
            end               
        end
    end
    
end
