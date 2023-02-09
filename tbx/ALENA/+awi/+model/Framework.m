classdef (ConstructOnLoad) Framework < awi.model.Entity ... is itself an entity belonging to the AWI Framework
        & mvc.model.Application ...                         %   AND it an Application in its own right
        & mvc.mixin.Analysable  ...                         %   supporting various types of analyses
        & mvc.mixin.Auditable
    %& mvc.mixin.Pathadjustable
    %Framework Class defining the top level object for the AWI Framework.
    %
    % An interactive environment for aero-elastic analysis, it is also a
    % component of that analysis in its own right (awi.model.Entity) and it
    % is also an Application (mvc.model.Application). Finally, it supports
    % the analysis of its contents (mvc.mixin.Analysable).
    %
    % Defines
    %   * Collection Specification:
    %       -
    %   * Import Specification:
    %       -
    %   * Export Specification:
    %       -
    %   * Context Menus:
    %       -
    
    properties (Dependent) % Aircraft, LoadCases, ResultsSets
        Aircraft;       % Handy shortcut
        LoadCases;      % Handy shortcut
        ResultSets;     % Handy shortcut
        ControlSurfaces % Handy shortcut
    end
    
    properties (Dependent)
        
        %Further shortcuts, to help with logic (mostly internal, could consider making these protected)
        HasAircraft;
        HasVisibleAircraft;
        HasControlSurfaces;
        
    end
    
    methods % set / get
        
        function val = get.ControlSurfaces(obj)
            if obj.HasAircraft
                val = obj.Aircraft.ControlSurfaces;
            else
                val = [];
            end
        end
        
        function val = get.HasControlSurfaces(obj)
            val = ~isempty(obj.ControlSurfaces);
        end
        
        function val = get.HasVisibleAircraft(obj)
            
            %Start here
            val = obj.HasAircraft;
            
            %Anything ?
            if ~val
                
                %No
                return;
                
            end
            
            %But does the Aircraft actually have any visible content ?
            val = ~isempty(findall(obj.Aircraft, 'isa', 'mvc.mixin.Drawable', '-and', 'Visible', true));
            
        end
        
        function val = get.HasAircraft(obj)
            
            %Simple
            val = ~isempty(obj.Aircraft);
            
        end
        
        function set.Aircraft(obj, val)
            
            %Ensure built before adding to framework
            if isa(val, 'awi.mixin.Buildable')
                build(val);
            end
            
            %Check for existing aircraft
            if isempty(obj.Children)
                
                %No need to check
                idx = [];
                
            else
                
                %Use find to locate
                idx = find(obj.Children, '-index', 'isa', 'awi.model.Aircraft');
                
            end
            
            %Anything ?
            if isempty(idx)
                
                %Allow for setting no aircraft
                if isempty(val)
                    
                    %Nothing to do
                    
                else
                    
                    %Just add it
                    obj.add(val);
                    
                end
                
            else
                
                %Allow for setting no aircraft
                if isempty(val)
                    
                    %Remove it
                    obj.Children(idx) = []; % .removeThis;
                    
                else
                    
                    %Replace it
                    obj.Children(idx) = val;
                    
                end
                
            end
            
        end
        
        function val = get.Aircraft(obj)
            
            %Anything ?
            if isempty(obj.Children)
                
                %No
                val = [];
                
            else
                
                %Use find to locate
                val = find(obj.Children, 'isa', 'awi.model.Aircraft');
                
            end
            
        end
        
        function set.LoadCases(obj, val)
            
            %Find collector (TODO: is this safe ?)
            col = findall(obj.Children, 'Name', 'Load Cases');
            
            %If nothing
            if isempty(col)
                
                %Anything to add ?
                if isempty(val)
                    
                    %No
                    
                else
                    
                    %Just add whatever
                    obj.add(val);
                    
                end
                
            else
                
                %Assign
                col.Children = val;
                
            end
            
        end
        
        function val = get.LoadCases(obj)
            
            %Anything ?
            if isempty(obj.Children)
                
                %No
                val = [];
                
            else
                
                %Use find to locate
                val = findall(obj.Children, 'isa', 'awi.model.LoadCase');
                
            end
            
        end
        
        function set.ResultSets(obj, val)
            
            %Find collector (TODO: is this safe ?)
            col = findall(obj.Children, 'Name', 'Result Sets');
            
            %If nothing
            if isempty(col)
                
                %Anything to add ?
                if isempty(val)
                    
                    %No
                    
                else
                    
                    %Just add whatever
                    obj.add(val);
                    
                end
                
            else
                
                %Assign
                col.Children = val;
                
            end
            
        end
        
        function val = get.ResultSets(obj)
            
            %Anything ?
            if isempty(obj.Children)
                
                %No
                val = [];
                
            else
                
                %Use find to locate
                val = findall(obj.Children, 'isa', 'mvc.model.ResultSet');
                
            end
            
        end
        
    end
    
    methods % construction / destruction
        
        function obj = Framework(varargin)
            
            %Pass it on
            obj@awi.model.Entity('Name', 'ALENA Framework', varargin{:});
            
            %Configure collectables - the Aircraft does not directly parent these
            obj.rmCollectionSpec('Bluff Body', 'Lifting Surface', 'Point Mass', 'Coordinate System');
            
            %But it does collect
            obj.addCollectionSpec(@awi.model.Aircraft, 'Aircraft', [], 1); %No more than one of these
            obj.addCollectionSpec(@awi.model.LoadCase, 'Load Case', [], [], [], [], 'Load Cases'); %Grouped
            obj.addCollectionSpec(@awi.model.BeamModel, [], [], [], true, [], 'Result Sets'); %Hidden, grouped
            %obj.addCollectionSpec(@awi.model.ResultSet, [], [], [], true, [], 'Result Sets'); %Hidden, grouped
%             obj.addCollectionSpec(@awi.model.Material, 'Material', [], [], [], [], 'Materials'); %Grouped
%             obj.addCollectionSpec(@awi.model.CrossSection, 'CrossSection', [], [], [], [], 'CrossSections'); %Grouped

            %Framework does not actually draw itself, just its content
            obj.Visible = false;
            
            %Framework looks best if we do NOT impose uniqueness on all children, regardless of class
            obj.AllNamesUnique = false;
            
            %Configure serialize spec
            obj.SerializeSpec = [{'*.awi', 'ALENA Framework files (*.awi)', [], []}; ...
                obj.SerializeSpec];
            
            %Configure importability, by tweaking the XML import behaviour
            obj.ImportXMLDefaultClass = 'awi.model.Entity';
            
            %And adding capability to import from FAME
            obj.ImportSpec(end+1:end+2,:) = ...
                { ...
                '*.fm4', 'FAME Input files (*.fm4)'     , @import_fm4; ...
                '*.awi', 'ALENA Framework files (*.awi)', @import_mat};
            
            %If debuggable
            if isa(obj, 'mvc.mixin.Debugable')
                
                %Include self-test as an option
                addContext(obj, '|Debug>Selftest...', 'selftest');
                
            end
            
            %Configure analyses
            obj.addAnalysis(@generateNastranModel , @canConvertToFE , 'Generate FE model>MSC.Nastran'); %-> Is this the best way to do this? Could we just export to these file formats instead? I think that would be better...
%             obj.addAnalysis(@generateNeoCassModel , @canConvertToFE , 'Generate FE model>NeoCass');
%             obj.addAnalysis(@generateAeroFlexModel, @canConvertToFE , 'Generate FE model>AeroFlex');
            obj.addAnalysis(@massAnalysis   , @canMassAnalysis   , 'Model Checking>Mass');
            %obj.addAnalysis(@massAnalysis, @canMassAnalysis, 'Model Checking>Element Connectivity');
            obj.addAnalysis(@trimAnalysis   , @canMassAnalysis   , 'Static Aeroelastic>Trimmed Manoeuvre'); %TODO - Remove/finesse this
            obj.addAnalysis(@gustAnalysis   , @canMassAnalysis   , 'Dynamic Aeroelastic>Dynamic Gust');  
            obj.addAnalysis(@flutterAnalysis, @canFlutterAnalysis, 'Dynamic Aeroelastic>Flutter (Future Work) ');                      
            obj.addAnalysis(@sizeAnalysis   , @canMassAnalysis   , 'Sizing');

            %If caller does not want it back
            if nargout == 0
                
                %Open in corresponding view
                awi.view.Framework(obj);
                
            end
            
        end
        
    end
    
    methods % Selftest
        
        function varargout = selftest(obj, varargin)
            
            %Pass it on to test class
            [varargout{1:nargout}] = awi.model.SelfTest(obj, varargin{:});
            
        end
        
    end
    
    methods % Analyses
        
        function b = canTrimAnalysis(obj)
            
            %If we have an aircraft, we can trim it,
            % but only if we have one or more load cases
            b = ~isempty(obj.Aircraft) && ~isempty(obj.LoadCases);
            
        end
        
        function b = canSizeAnalysis(obj)
            
            %If we have an aircraft, we can trim it,
            % but only if we have one or more load cases
            b = ~isempty(obj.Aircraft) && ~isempty(obj.LoadCases);
            
        end
        
        function b = canFlutterAnalysis(obj)
            
            b = false;
            return
            
            %If we have an aircraft, we can do flutter analysis
            b = ~isempty(obj.Aircraft);
            
            %Or can we ?  Not sure.  For time being just say no
            b = false;
            
        end
        
        function b = canMassAnalysis(obj)
            
            %If we have an aircraft, we can analyse its mass
            b = ~isempty(obj.Aircraft);
            
        end
        
        function b = canConvertToFE(obj)
            
            %If we have an aircraft, we can convert it to a FE model as it
            %will inherit from 'awi.mixin.FEable'
            b = ~isempty(obj.Aircraft);
        end
        
        function varargout = massAnalysis(obj, varargin)
            %Doing Mass Analysis is same as creating a new view of type Mass Analysis
            
            %Look for a view manager associated with this object
            mgr = viewManager(obj, varargin{:});
            
            %Anything ?
            if isempty(mgr)
                
                %No, create a stand-alone view
                [varargout{1:nargout}] = awi.view.MassAnalysis(obj);
                
            else
                
                %Pass it on
                [varargout{1:nargout}] = mgr.addView(@awi.view.MassAnalysis);
                
            end
            
        end
        
        function varargout = trimAnalysis(obj, varargin)
            %Doing Trim Analysis is same as creating a new view of type Trim Analysis
            
            %Look for a view manager associated with this object
            mgr = viewManager(obj, varargin{:});
            
            %Anything ?
            if isempty(mgr)
                
                %No, create a stand-alone view
                [varargout{1:nargout}] = awi.view.TrimAnalysis(obj);
                
            else
                
                %Pass it on
                [varargout{1:nargout}] = mgr.addView(@awi.view.TrimAnalysis);
                
            end
            
        end
        
        function varargout = gustAnalysis(obj, varargin)
            %Doing Loac Case Analysis is same as creating a new view of type Trim Analysis
            
            %Look for a view manager associated with this object
            mgr = viewManager(obj, varargin{:});
            
            %Anything ?
            if isempty(mgr)
                
                %No, create a stand-alone view
                [varargout{1:nargout}] = awi.view.GustAnalysis(obj);
                
            else
                
                %Pass it on
                [varargout{1:nargout}] = mgr.addView(@awi.view.GustAnalysis);
                
            end
            
        end
        
        function varargout = sizeAnalysis(obj, varargin)
            %Doing Size Analysis is same as creating a new view of type Size Analysis
            
            %Look for a view manager associated with this object
            mgr = viewManager(obj, varargin{:});
            
            %Anything ?
            if isempty(mgr)
                
                %No, create a stand-alone view
                [varargout{1:nargout}] = awi.view.SizeAnalysis(obj);
                
            else
                
                %Pass it on
                [varargout{1:nargout}] = mgr.addView(@awi.view.SizeAnalysis);
                
            end
            
        end
        
        function generateNastranModel(obj)
            
            %Do we already have an instance of MSC.Nastran FEM in the
            %collection?
            mscFEM = findall(obj, 'Name', 'MSC.Nastran FEM');
            
            if isempty(mscFEM)
                name = 'MSC.Nastran FEM';
            else
                nMSC = numel(mscFEM);
                name = ['MSC.Nastran FEM_#', num2str(nMSC + 1)];
            end
            
            %Make a generic beam model (data structure)
            BeamModel = convertToFE(obj.Aircraft);

            %Insert an intermediate BeamModel
            B = awi.model.BeamModel( ...
                'Aircraft', obj.Aircraft.Name, ...
                'Name', name, ...
                'BM', BeamModel);
            
            %Add it to the model
            obj.add(B);
        end
        
        function generateNeoCassModel(obj)
            
            %Do we already have an instance of NeoCass FEM in the
            %collection?
            neoFEM = findall(obj, 'Name', 'NeoCass FEM');
            
            if isempty(neoFEM)
                name = 'NeoCass FEM';
            else
                nNeo = numel(neoFEM);
                name = ['NeoCass FEM_#', num2str(nNeo + 1)];
            end
            
            %Make a basic beam model (structure)
            BeamModel = convertToFE(obj.Aircraft);
            
            %Convert to a specific NeoCass model
            
            
            %Insert an intermediate Beam Model
            B = awi.model.BeamModel( ...
                'Aircraft', obj.Aircraft.Name, ...
                'Name', name);
            
            %Add it to the model
            obj.add(B);
        end
        
        function generateAeroFlexModel(obj)
            
        end
        
    end
    
    methods % Custom importers
        
        function newobj = import(obj, varargin)
            
            %Look for a "force" flag in varargin
            bForce = false;
            idx = find(cellfun(@(x)ischar(x) && strcmpi('force', x), varargin));
            if ~isempty(idx)
                bForce = varargin{idx+1};
                varargin(idx + [0,1]) = [];
            end
            
            %Start with base-class
            [newobj, logfcn] = import@mvc.mixin.Importable(obj, varargin{:});
                                   
            %Cancelled ?
            if isempty(newobj)
                
                %Done with progress window
                if ~isempty(logfcn)
                    logfcn('close');
                end
                
                %Bail
                return;
                
            elseif nargout
                
                %Done with progress window
                logfcn('close');
                
                %Caller wants the imported content back, so do nothing more here
                return;
                
            end
            
            if ~isempty(newobj) && strcmpi(newobj.Name, 'CPACS') %CPACS?
                cpacs = awi.model.cpacs;
                cpacs.add(newobj.Children);
                newobj = convertCPACSToAWI(cpacs);
            end
            
            %Reformat hierarchy to store Material/Profile data in the
            %correct place.
            ch   = newobj.Children;
            bMat = arrayfun(@(o) isa(o, 'awi.model.Material'), ch);
            Materials = findall(ch(~bMat), 'Type', 'Material');
            if ~isempty(Materials)
                Materials = detach(Materials);
                add(newobj, Materials);
            end
            
            %For audit trail
            msg = {};
            
            %What have we got ?  Check whether we've just imported an entire framework session,
            if isa(newobj, class(obj)) ...                              % An entire Framework
                    || isa(newobj, 'awi.model.Aircraft') ...            % OR just an Aircraft
                    || isa([newobj.Children], 'awi.model.LoadCase') ... % OR just a collection of load cases
                    || isa(newobj, 'awi.model.LoadCase') ...            % OR just one load case (or a load case array)
                    || isa([newobj.Children], 'awi.model.ResultSet') ...% OR just a collection of result sets
                    || isa(newobj, 'awi.model.ResultSet')               % OR just one result set (or a result set array)
                
                %Deal with the Aircraft first
                AC = findall(newobj, 'isa', 'awi.model.Aircraft');
                if ~isempty(AC)
                    
                    %Not impossible
                    assert(numel(AC) == 1, 'framework permits only a single Aircraft model in each session');
                    
                    %If no description
                    if isempty(AC.Description)
                        
                        %Add something suitable (TODO: use Auditable if possible, rather than description ?)
                        AC.Description = ['Aircraft imported from file ''', obj.ImportFullFile, ''''];
                        
                    end
                    
                    %If force flag is set by caller, or if we have no extant aircraft, or if user confirms
                    if bForce || isempty(obj.Aircraft) || obj.confirm('Overwrite existing Aircraft ?')
                        
                        %Assign in session
                        obj.Aircraft = AC;
                        
                        %Initialise any imported Aircraft to be structurally LOCKED
                        set(obj.Aircraft, 'StructureLocked', true);
                        
                        %For audit trail
                        msg{end+1} = 'aircraft';
                        
                    else
                        
                        %Need to throw the imported aircraft away, else we might pick it up again later
                        AC.removeThis;
                        
                    end
                    
                end
                
                %If we have one or more Load Cases
                LC = findall(newobj, 'isa', 'awi.model.LoadCase');
                if ~isempty(LC)
                    
                    %If no description
                    b = cellfun(@isempty, {LC.Description});
                    [LC(b).Description] = deal(['Loadcase imported from file ''', obj.ImportFullFile, '''']);
                    
                    %Overwrite any existing load-cases, or append ?
                    %If force flag is set by caller, or if we have no extant load-cases
                    if bForce || isempty(obj.LoadCases)
                        
                        %Just go with overwrite
                        obj.LoadCases = LC;
                        
                        %For audit trail
                        if numel(LC) == 1
                            msg{end+1} = 'loadcase';
                        else
                            msg{end+1} = [num2str(numel(LC)), ' loadcases'];
                        end
                        
                    else
                        
                        %Ask the user
                        switch choose(obj, 'Overwrite extant Load Case(s), or append with those imported from file ?', ...
                                'Overwrite', 'Append', 'Do neither', 'Overwrite')
                            
                            case 'Overwrite'
                                
                                %Do it
                                obj.LoadCases = LC;
                                
                                %For audit trail
                                if numel(LC) == 1
                                    msg{end+1} = 'loadcase';
                                else
                                    msg{end+1} = [num2str(numel(LC)), ' loadcases'];
                                end
                                
                            case 'Append'
                                
                                %Do it
                                obj.add(LC);
                                
                                %For audit trail
                                if numel(LC) == 1
                                    msg{end+1} = 'loadcase';
                                else
                                    msg{end+1} = [num2str(numel(LC)), ' loadcases'];
                                end
                                
                            otherwise
                                
                                %Need to throw the imported loadcases away, else we might pick it up again later
                                LC.removeThis;
                                
                        end
                        
                    end
                    
                end
                
                %If we have one or more Result Sets
                %                 RS = findall(newobj, 'isa', 'awi.model.ResultSet');
                RS = findall(newobj, 'isa', 'mvc.model.Collector', 'Name', 'Result Sets');
                
                if ~isempty(RS)
                    
                    %                     %If no description
                    %                     b = cellfun(@isempty, {RS.Description});
                    %                     [RS(b).Description] = deal(['Resultset imported from file ''', obj.ImportFullFile, '''']);
                    %
                    %Overwrite any existing result sets or append ?
                    %If force flag is set by caller, or if we have no extant result sets
                    if bForce || isempty(obj.ResultSets)
                        
                        %Just go with overwrite
                        obj.ResultSets = RS.Children;
                        
                        %For audit trail
                        if numel(RS.Children) == 1
                            msg{end+1} = 'result';
                        else
                            msg{end+1} = [num2str(numel(RS.Children)), ' results'];
                        end
                        
                    else
                        
                        %Ask the user
                        switch choose(obj, 'Overwrite extant Result Sets(s), or append with those imported from file ?', ...
                                'Overwrite', 'Append', 'Do neither', 'Overwrite')
                            
                            case 'Overwrite'
                                
                                %Do it
                                obj.ResultSets = RS.Children;
                                
                                %For audit trail
                                if numel(RS.Children) == 1
                                    msg{end+1} = 'result';
                                else
                                    msg{end+1} = [num2str(numel(RS.Children)), ' results'];
                                end
                                
                            case 'Append'
                                
                                %Do it
                                obj.add(RS.Children);
                                
                                %For audit trail
                                if numel(RS.Children) == 1
                                    msg{end+1} = 'result (appended to existing)';
                                else
                                    msg{end+1} = [num2str(numel(RS.Children)), ' results (appended to existing)'];
                                end
                                
                            otherwise
                                
                                %Need to throw the imported result-sets away, else we might pick it up again later
                                RS.Children.removeThis;
                                
                        end
                        
                    end
                    
                end
                
                %Eliminate any collectors that have no content
                newobj.Children(arrayfun(@(x)isa(x, 'mvc.model.Collector') && isempty(x.Children), newobj.Children)) = [];
                
                %And if newobj has any immediate children that are none of the above
                if ~isempty(newobj.Children)
                    
                    %Add them at top-level
                    obj.add(newobj.Children);
                    
                end
                
            elseif isa(newobj, 'mvc.mixin.Collectable')
                
                %To where could the item(s) be moved ?
                lst = obj.Root.flatlist(2);
                b = arrayfun(@(x)x.canAdd(newobj), lst);
                assert(any(b), 'no feasible destinations for imported content');
                lst = lst(b);
                
                %Only one ?
                if numel(lst) == 1
                    
                    %Go for it
                    dest = lst;
                    
                else
                    
                    %User selects destination
                    dest = obj.selectFrom(lst, 'Add imported content to...', ...
                        'SelectionMode', 'single');
                    
                    %Cancelled ?
                    if isempty(dest)
                        
                        %Bail out
                        newobj = [];
                        return;
                        
                    end
                    
                end
                
                %Do it
                dest.add(newobj);
                
                %For audit trail
                if numel(newobj) == 1
                    msg{end+1} = class(newobj);
                else
                    msg(end+1:end+numel(newobj)) = cellfun(@class, newobj, 'UniformOutput', false);
                end
                
            end
            
            %Did we actually import anything ?
            if isempty(msg)
                
                %No
                
            elseif isa(obj, 'mvc.mixin.Auditable')
                
                %Update audit trail, if appropriate, taking details from helper object
                if isa(newobj, 'mvc.mixin.Auditable')
                    addAuditTrailEntry(obj, newobj.AuditTrail);
                end
                
                %And summarise what actually happened to this object
                addAuditTrailEntry(obj, ['Imported ', strjoin(msg, ', '), ' from ''', obj.ImportFullFile, '''']);
                
            end
            
            %Done with progress window
            logfcn('close');
            
        end
        
        function newobj = import_fm4(obj, logfcn, varargin)
            
            %Implement ui as standard
            isBatch = false;
            
            %Determine if FAME import is being called in 'batch' mode.
            if any(contains(varargin, '-batch'))
                isBatch = true;
            end
            
            %Find the folder containing all FAME files
            fameFolder = i_findFameFolder(obj.ImportFullFile);
            
            function fameFolder = i_findFameFolder(importFile)
                %i_findFameFolder Returns the parent directory of the FAME
                %files.
                
                %Identifier for finding the correct path
                bWeightsFound = false;
                
                %Initialise the pathname as dn
                dn = fileparts(importFile);
                
                %Loop until its found
                while ~bWeightsFound
                    dn_ = fileparts(dn);
                    if ~isempty(regexp(dn_, '06_weights$', 'ONCE'))
                        bWeightsFound = true;
                    end
                    if strcmp(dn,dn_)
                        %logfcn('**WARNING** 06_weights folder not found');
                        break;
                    end
                    dn = dn_;
                end
                
                fameFolder = fileparts(dn);
                
            end
            
            %% Ask the user what files they want to import 
            
            % Create the Fame Model object
            a_fame = Fame2mat;
            
            %Only run the import manager if we are not in '-batch' mode
            if ~isBatch
                
                %Configure the import files
                importFileOpts = i_configureImportFiles;
                
                %Check for special tokens
                if ischar(importFileOpts) && strcmp(importFileOpts, '-QuitImport')
                    newobj = [];
                    logfcn('User terminated the import process');
                    return
                end
                
                %Assign to the 'fame2mat' object
                set(a_fame, importFileOpts(1, :), importFileOpts(2, :));
                
            end
            
            function pav = i_configureImportFiles
                %i_configureImportFiles Opens up a window that asks the
                %user to detail which of the FAME input files they have
                %access to. 
                %
                % This allows the user to skip parts of the of the import
                % process that would usually ask for a file that has not
                % been automatically found.
                
                %FAME files/directories that the user can select
                fameFiles = { ...
                    'Fame output directory (.pv8)', 'HasPV8'    ; ...
                    'Aircraft geometry (.tux_xml)', 'HasTuxXML' ; ...
                    'Mass distribution (.205)'    , 'Has205'    ; ...
                    'Fuel mass files'             , 'HasFuel'  };
                
                %Defaults for the selection window
                spc = 10;  %[px], spacing around elements in the layout
                pad = 5;   %[px], padding between elements in each layout
                w   = 500; %[px], width of the window
                h   = 200; %[px], height of the window
                bh  = 20;  %[px], height of the uicontrol objects
                th  = 50;  %[px], height of the text uicontrol
                
                %Position the window in the middle of the screen
                monPos = get(0, 'MonitorPositions');
                monPos = monPos(1, :);                
                winPos = [ ...
                    monPos(1) + monPos(3) / 2 - w / 2, ...
                    monPos(2) + monPos(4) / 2 - h / 2, ...
                    w, h];
                
                %Make the graphics objects
                hSel = figure('Name', 'Configure FAME import files', ...
                    'NumberTitle', 'off'   , ...
                    'MenuBar'    , 'none'  , ...
                    'ToolBar'    , 'none'  , ...
                    'Units'      , 'pixels', ...
                    'Position'   , winPos);
                vBox     = uix.VBox('Parent', hSel, 'Spacing', spc);
                vBBox(1) = uix.VButtonBox('Parent', vBox, ...
                    'Spacing', spc, 'Padding', pad, ...
                    'HorizontalAlignment', 'left');
                hBox     = uix.HBox('Parent', vBox, 'Spacing', spc);
                hScroll  = uix.ScrollingPanel('Parent', hBox);
                vBBox(2) = uix.VButtonBox('Parent', hScroll, 'Spacing', spc, ...
                    'HorizontalAlignment', 'left', 'Padding', pad);
                vBBox(3) = uix.VButtonBox('Parent', hBox, 'Spacing', spc, ...
                    'Padding', pad);
                vBox.Heights = [th -1];
                hBox.Widths  = [-4 -1];
                uicontrol('Parent', vBBox(1) , 'Style', 'text' , ...
                    'String', 'Which FAME files/directories do you have access to?', ...
                    'HorizontalAlignment', 'left');
                
                %Push buttons for progressing the code
                hImp = uicontrol('Parent', vBBox(3), 'Style', 'push', ...
                    'String', 'Import!', 'Callback', @runImport, 'Enable', 'off');
                uicontrol('Parent', vBBox(3), 'Style', 'push', ...
                    'String', 'Cancel', 'Callback', @cancelImport);
                
                %Make the check boxes
                hOpt = gobjects(1, size(fameFiles, 1) + 1);
                for iOpt = 1 : size(fameFiles, 1)
                    hOpt(iOpt) = uicontrol('Parent' , vBBox(2), ...
                        'Style'   , 'check'           , ...
                        'String'  , fameFiles{iOpt, 1}, ...
                        'Tag'     , fameFiles{iOpt, 2});
                end
                
                %Make additional option for if the user doesn't know
                hOpt(5)  = uicontrol('Parent', vBBox(2), 'Style', 'check', ...
                    'String', 'Not sure', 'Callback', {@clearOpts, hOpt(1 : 4), hImp});
                set(hOpt(1 : end - 1), 'Callback', {@clearOpts, hOpt(end), hImp})   
                
                %Update the layout dimensions
                hRatio = abs(hBox.Widths) ./ sum(abs(hBox.Widths));
                vBBox(1).ButtonSize = [0.9 * w, th];
                vBBox(2).ButtonSize = [(hRatio(1) - 0.1) * w, bh];
                hScroll.MinimumHeights = numel(vBBox(2).Children) * vBBox(2).ButtonSize(2) + ...
                    vBBox(2).Spacing + (numel(vBBox(2).Children) - 1) * pad;
                
                %Wait for user input
                uiwait(hSel);
                                
                function clearOpts(~, ~, hOpt, hImp)
                    %clearOpts Sets all uicontrol objects in 'hOpt' to have
                    %'Value' = 0 and enables the import button 'hImp'.
                    
                    %Uncheck all options
                    set(hOpt, 'Value', 0);
                    
                    %Enable the import
                    hImp.Enable = 'on';
                    
                end
                
                function cancelImport(src, ~)
                    %cancelImport Tags the figure so the import method
                    %knows to quit out of the import process.
                    
                    hFig = ancestor(src, 'Figure');
                    hFig.Tag = 'QuitImport';
                    uiresume(hFig);
                    
                end
                
                function runImport(src, ~)
                    %runImport Sets 'uiresume' on the current figure.
                    
                    hFig = ancestor(src, 'Figure');
                    uiresume(hFig);
                    
                end
                
                %Check to see if the user has close the window already
                if ~isvalid(hSel)
                    pav = '-QuitImport';
                    return
                end
                
                %Has the figure been marked?
                if ~isempty(hSel.Tag)
                    switch hSel.Tag
                        case 'QuitImport'
                            pav = '-QuitImport';  
                            close(hSel)
                            return
                        otherwise
                    end
                end
                
                %Get import options
                prop  = {hOpt.Tag};
                val   = {hOpt.Value};
                index = cellfun(@isempty, prop);
                pav   = [prop(~index) ; val(~index)];
                
                %If the user isn't sure then flag everything as available
                %so the FAME import object will allow the user to manually
                %select them.
                if hOpt(end).Value
                   [pav(2, :)] = deal({true});  
                end
                
                %Close the figure
                close(hSel);
                
            end
                        
            %% Import Beam Model if found
            
            logfcn('Generating an Aeroelastic Beam Model of the FAME model: ');            

            % Assign the input file to the Fame object:
            a_fame.Inp.fm4File{1} = obj.ImportFullFile;
            
            % Redefine the directory used by the Fame Tool
            a_fame.Dirs.MainFolder = fameFolder;
            
            a_fame = run(a_fame, logfcn);
            
            % Recall the fm4Input
            logfcn('Converting FAME file to AWI Framework ...');
            [Aircraft_, FuelCases] = FameModel2AWIFramework(a_fame.Fame, a_fame.Mdl, logfcn);
            logfcn('Conversion Complete.');
                        
            logfcn('Building the AWI Framework model ...');
            %Attempt to build the model using any 'Buildable' properties
            idx = arrayfun(@(x) isa(x, 'awi.mixin.Buildable'), Aircraft_.Children);
            arrayfun(@build, Aircraft_.Children(idx));
            logfcn('Build Complete!');
            
            %Create a new instance of 'awi.model.Framework' and add 
            %'Aircraft_' as a child of this new instance.
            newobj = obj.new('Children', Aircraft_);
            
            %% Import Load Cases if found
            
            if ~isempty(a_fame.Fame.LoadCases)
                
                logfcn(['Extracting the Fame Load Cases: ' a_fame.Inp.loadCaseFiles{1}]);
                
                %Index the load cases from the FAME structure & the .fm4 
                lc = a_fame.Fame.LoadCases;
                lc_fm4 = [];
                if isfield(a_fame.Fame.fm4Input, 'WING') && ...
                        isfield(a_fame.Fame.fm4Input.WING, 'LOAD_CASES')
                    lc_fm4 = a_fame.Fame.fm4Input.WING.LOAD_CASES;                    
                end
                    
                logfcn('Converting the Fame Load Cases to an AWI framework load case...');
                %Convert the load case structures into 'awi.model.LoadCase'
                %objects
                LC = awi.model.LoadCase.fame2obj(lc, lc_fm4, FuelCases);
                logfcn('Load case conversion complete!');
                
                logfcn('Assigning control surface deflections ...');
                %Grab the control surface names
                %   - Because some of the control surfaces are associated
                %   with collector nodes we will need to flatlist the
                %   collection to find all control surfaces...
                collection = flatlist(newobj);
                cidx   = arrayfun(@(o) isa(o, 'awi.model.ControlSurface'), collection);
                cs     = collection(cidx);
                csName = {cs.Name};
                
                % Find the spoilers
                %   - TODO : This needs to be improved! Could I use the new
                %   HingePosition property of the control surface?
                                
                % splIdx = ~cellfun(@isempty,regexp(csName,'spo'));
                % ailIdx = ~cellfun(@isempty,regexp(csName,'ail'));
                
                % Get the starboard control surfaces
                splIdx  = strncmpi('spo',csName,3);
                ailIdx  = strncmpi('ail',csName,3);
                elevIdx = strncmpi('ele',csName,3);
                
                % Get the port side control surfaes
                port_splIdx  = strncmpi('pspo',csName,4);
                port_ailIdx  = strncmpi('pail',csName,4);
                port_elevIdx = strncmpi('pele',csName,4);
                
                %Index for all spoilers/ailerons on the port and starboard
                %wing
                csFixedIdx = or(or(splIdx, ailIdx), or(port_splIdx, port_ailIdx));
                csFreeIdx = or(elevIdx, port_elevIdx);
                
                AileronFactor = 1.5; % I don't know why but this seems to be needed!
                
                %Assign control surface deflections to the
                for i = 1:numel(lc)
                    %Assign values
                    LC(i).CsDeflection(splIdx) = 0.5*lc(i).spoiler;
                    LC(i).CsDeflection(ailIdx) = AileronFactor*lc(i).aileron;
                    LC(i).CsDeflection(port_splIdx) = 0.5*lc(i).spoiler;
                    LC(i).CsDeflection(port_ailIdx) = AileronFactor*lc(i).aileron;
                    %Set control surface deflection type to 'fixed'
                    LC(i).CsDeflecType(csFixedIdx) = {'fixed'};
                    LC(i).CsDeflecType(csFreeIdx) = {'free'};
                    
                end
                logfcn('Control surface deflections assigned!')
                
                logfcn('Adding the load cases to the collection ...')
                %Add the load case objects to the model
                newobj.add(LC);
                logfcn('Load cases added to the collection!')
                
                logfcn('Load case extraction complete!');
                
            end
            
            %% Convert the FAME beam model to an AWI beam model
            
            if numel(a_fame.Mdl.Grid)
                
                logfcn('Converting the fame object to a structure ...');
                BeamModel = FameObj2BeamModel(a_fame.Mdl);
                logfcn('Beam model structure generated!');
                
                %Insert an intermediate Beam Model
                B = awi.model.BeamModel( ...
                    'Aircraft', newobj.Aircraft.Name, ...
                    'LoadCase', {newobj.LoadCases.Name},...
                    'BM'      , BeamModel,...
                    'Name'    ,'FAME FEM');
                
                %Create one 'AeroelasticResult' for each load case in the FAME
                %result set
                nResults = numel(a_fame.Fame.Results.InternalLoads);
                StbdWingResults = arrayfun(@(i) awi.model.AeroelasticResult, 1 : nResults, 'Unif', false);
                StbdWingResults = horzcat(StbdWingResults{:});
                TrimResults = arrayfun(@(i) awi.model.ResultSet, 1 : nResults, 'Unif', false);
                TrimResults = horzcat(TrimResults{:});
                
                %Only one sizing result is required per FAME file
                FameSizeResult = awi.model.SizingResult;
                SizeResult     = awi.model.ResultSet('Name', 'FAME Sizing Result');
                
                %Set names of the 'TrimResults' to match meaningful names from
                %FAME input
                LCi = [a_fame.Fame.LoadCases.LC]; %Loadcase number
                Labels = arrayfun(@(i) sprintf('FAME Result (LC %i)', i), LCi, 'Unif', false);
                set(TrimResults, {'Name'}, Labels');
                
                %FAME  results are only for the starboard wing. Find this
                %object in the collection and pass its handle to 'TrimResults'.
                %   - Assume the name of the starboard wing object is
                %   'StbdWing'. 'findall' automatically returns a default
                %   empty object and if we pass this into 'TrimResults' it will
                %   be thrown away by the method 'set.Beam'.
                StbdWing = findall(Aircraft_, 'Name', 'StbdWing');
                assert(~isempty(StbdWing), 'No object called ''StbdWing'' found in the AWI model');
                set(StbdWingResults, 'Beam', StbdWing);
                set(FameSizeResult , 'Beam', StbdWing);
                
                %Format the results from FAME
                %   - The eta position is the same for all results of the
                %   same type. i.e. The eta position for all internal loads
                %   is the same, however the aerodynamic results might have
                %   their own eta positions.
                %   - Displacements
                eta_disp = {a_fame.Fame.Results.Deformation.eta}';
                Tx = {a_fame.Fame.Results.Deformation.Ux}';
                Tz = {a_fame.Fame.Results.Deformation.Uz}';
                %   - Internal Loads
                eta_loads = {a_fame.Fame.Results.InternalLoads.eta}';
                Fx  = {a_fame.Fame.Results.InternalLoads.Fx}';
                Fy  = {a_fame.Fame.Results.InternalLoads.Fy}';
                Fz  = {a_fame.Fame.Results.InternalLoads.Fz}';
                Mx  = {a_fame.Fame.Results.InternalLoads.Mx}';
                My  = {a_fame.Fame.Results.InternalLoads.My}';
                Mz  = {a_fame.Fame.Results.InternalLoads.Mz}';
                %   - Aerodynamic Results
                eta_aero = {a_fame.Fame.Results.Aerodynamic.eta}';
                Cl       = {a_fame.Fame.Results.Aerodynamic.Cl}';
                Cl_l     = {a_fame.Fame.Results.Aerodynamic.Cl_l}';
                Cm       = {a_fame.Fame.Results.Aerodynamic.Cm}';
                
                %Assign to the objects
                set(StbdWingResults, {'FxFAME_eta'}, eta_loads);
                set(StbdWingResults, {'FxFAME'}    , Fx);
                set(StbdWingResults, {'FyFAME_eta'}, eta_loads);
                set(StbdWingResults, {'FyFAME'}    , Fy);
                set(StbdWingResults, {'FzFAME_eta'}, eta_loads);
                set(StbdWingResults, {'FzFAME'}    , Fz);
                set(StbdWingResults, {'MxFAME_eta'}, eta_loads);
                set(StbdWingResults, {'MxFAME'}    , Mx);
                set(StbdWingResults, {'MyFAME_eta'}, eta_loads);
                set(StbdWingResults, {'MyFAME'}    , My);
                set(StbdWingResults, {'MzFAME_eta'}, eta_loads);
                set(StbdWingResults, {'MzFAME'}    , Mz);
                %
                set(StbdWingResults, {'Tx'}    , Tx);
                set(StbdWingResults, {'Tx_eta'}, eta_disp);
                set(StbdWingResults, {'Tz'}    , Tz);
                set(StbdWingResults, {'Tz_eta'}, eta_disp);
                %
                set(StbdWingResults, {'Cl'}      , Cl);
                set(StbdWingResults, {'Cl_eta'}  , eta_aero);
                set(StbdWingResults, {'Cl_l'}    , Cl_l);
                set(StbdWingResults, {'Cl_l_eta'}, eta_aero);
                set(StbdWingResults, {'Cm'}      , Cm);
                set(StbdWingResults, {'Cm_eta'}  , eta_aero);
                
                %Aggregate the beam results into the 'ResultsSet' objects
                BR = num2cell(StbdWingResults);
                set(TrimResults, {'BeamResults'}, BR');
                
                %Assign the FAME sizing results data to the 'SizingResults'
                %object.
                if ~isempty(a_fame.Fame.Results.Thickness)
                    FameSizeResult.USknT       = a_fame.Fame.Results.Thickness.Uskn;
                    FameSizeResult.USknT_eta   = a_fame.Fame.Results.Thickness.eta;
                    FameSizeResult.LSknT       = a_fame.Fame.Results.Thickness.Lskn;
                    FameSizeResult.LSknT_eta   = a_fame.Fame.Results.Thickness.eta;
                    FameSizeResult.UShellT     = a_fame.Fame.Results.Thickness.Ush;
                    FameSizeResult.UShellT_eta = a_fame.Fame.Results.Thickness.eta;
                    FameSizeResult.LShellT     = a_fame.Fame.Results.Thickness.Lsh;
                    FameSizeResult.LShellT_eta = a_fame.Fame.Results.Thickness.eta;
                    FameSizeResult.FSparT      = a_fame.Fame.Results.Thickness.Fsp;
                    FameSizeResult.FSparT_eta  = a_fame.Fame.Results.Thickness.eta;
                    FameSizeResult.RSparT      = a_fame.Fame.Results.Thickness.Rsp;
                    FameSizeResult.RSparT_eta  = a_fame.Fame.Results.Thickness.eta;
                end
                
                %Aggregate the sizing results into the 'SizeResult' object
                SizeResult.BeamResults = FameSizeResult;
                
                %Add the trim  results and sizing results to the beam model
                B.add(TrimResults);
                B.add(SizeResult);
                
                %Add the 'awi.model.BeamModel' to the session
                newobj.add(B);
                
            end
            
            function fn = i_findfile(dn,mask)
                
                fn = [];
                
                d = dir(fullfile(dn,mask));
                if ~isempty(d)
                    
                    fn = fullfile(dn,d(1).name);
                    return;
                    
                else
                    d = dir(dn);
                    d(~[d.isdir]) = [];
                    d(ismember({d.name},{'.','..'})) = [];
                    for j = 1:numel(d)
                        fn = i_findfile(fullfile(dn,d(j).name),mask);
                        
                        if ~isempty(fn)
                            return
                        end
                    end
                end
                
            end
            
        end        
        
        function newobj = import_mat(obj, logfcn)
            
            newobj = [];
            
            return
            
            error('Import process for .mat files needs to be completed');
            
            %Treat as a generic .mat file
            temp = load(obj.ImportFullFile, '-mat');            
            
            %Parse the fields of 'temp' and decide what to do with the data
            
            %Action as appropriate
            
        end
        
    end
        
    methods (Static)% default objects etc.
        function Frmwrk = defaultFramework
            %defaultFramework Returns a default instance of the ALENA
            %Framework.
            %
            % Automatically populates the framework with a set of valid
            % objects, e.g. Aircraft, materials, etc.
            
            Frmwrk = awi.model.Framework('Name', 'ALENA Framework (default)');            
            AC = awi.model.Aircraft.defaultAircraft;
            Frmwrk.add(AC);
            
        end
    end
    
end

