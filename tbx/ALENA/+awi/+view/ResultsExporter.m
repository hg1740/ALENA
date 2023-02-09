classdef ResultsExporter < matlab.mixin.SetGet
    %ResultsExporter Handles the exporting of results from the AWI
    %Framework.
    
    %Primary properties
    properties
        %Handle to the view containing the results
        ResultsView
        %Selection in the results list box
        Selection = 1;
    end
        
    %Handles to graphics objects
    properties (SetAccess = protected)
        %Handle to the figure
        hFig        
        %Handle to the list of results
        hList        
        %Handle to the uicontrol object containing the path to the
        %user-specified results directory
        hResDir
        %Handle to the uicontrol object containing the list of beam names
        hBeamName
        %Handle to the uicontrol object containing the list of results
        %types
        hResType        
        %Handle to the checkbox for toggling the export to AWI file-format
        hAWI
        %Handle to the checkbox for toggling the export to CSV file-format
        hCSV
        %Handle to the checkbox for toggling the export to Excel file-format
        hExcel
    end
    
    %Shortcuts
    properties (Dependent)
        %Name of the directory wheree the results will be exported to
        ResDir
        %Results that are currently selected
        CurrentResults
        %Handle to all the 'BeamResults' objecets in the currently selected
        %results
        AllBeamResults
        %Handle to all the beam objects in the currently selected results
        AllBeams
        %Names of all the beams in the currently selected results
        AllBeamNames
        %Common beam names across all selected results
        CommonBeamNames
        %Names of all the common results types across the selection
        CommonResultTypes
        %Logical flag indicating whether the results should be exported in
        %the AWI file format
        ExportAWI
        %Logical flag indicating whether the results should be exported in
        %the CSV file format
        ExportCSV
        %Logical flag indicating whether the results should be exported in
        %the Excel file format
        ExportExcel
        %Combination of all export flags
        ExportFlag
        %Image data for the "File Open" icon (Output of 'imread')
        FileOpenIcon
    end
    
    %Helper properties for controlling the GUI appearance
    properties
        %Button height in pixels
        ButtonHeight = 26;
        %Number of pixels of extra space around the outside each GLT object
        Padding = 5;
        %Number of pixels of extra space to leave between each GLT object
        Spacing = 5;       
    end
    
    %Names of icon files
    properties (Constant)
        %Name of the "File Open" icon
        FileOpenIconName = 'file_open.png';
    end
    
    methods % set / get
        function set.ResultsView(obj, val)        %set.ResultsView
            %set.ResultsView Set method for the property 'ResultsView'
            %
            % - Must be an instance of 'awi.view.BeamResultsViewer'
            
            validateattributes(val, {'awi.view.BeamResultsViewer'}, ...
                {'scalar'}, class(obj), 'ResultsView');
            obj.ResultsView = val;
            
        end  
        function set.Selection(obj, val)          %set.Selection
            %set.Selection Set method for the property 'Selection'.
            
            %Sometimes an empty value can be passed through from the
            %ResultsViewer if it has just been initialised.
            if isempty(val) || ~isnumeric(val)
                val = 1;
            end
            
            validateattributes(val, {'numeric'}, {'integer', 'row', ...
                'positive'},  class(obj), 'Selection');
            obj.Selection = val;
            
        end
        function set.ButtonHeight(obj, val)       %set.ButtonHeight
           %set.ButtonHeight Set method for the property 'ButtonHeight'    
           
           validateattributes(val, {'numeric'}, {'scalar', 'integer', ...
               '>', 26}, class(obj), 'ButtonHeight');
           obj.ButtonHeight = val;
           
        end
        function set.Padding(obj, val)            %set.Padding
           %set.Padding Set method for the property 'Padding'    
           
           validateattributes(val, {'numeric'}, {'scalar', 'integer'}, ...
               class(obj), 'Padding');
           obj.Padding= val;
           
        end
        function set.Spacing(obj, val)            %set.Spacing
           %set.Spacing Set method for the property 'Spacing'    
           
           validateattributes(val, {'numeric'}, {'scalar', 'integer'}, ...
               class(obj), 'Spacing');
           obj.Spacing = val;
           
        end        
        function val = get.ResDir(obj)            %get.ResDir
            %get.ResDir Get method for the dependet property 'ResDir'.
            
            %Simple!
            val = obj.hResDir.String;
            
            %If it hasn't been defined then pass this back to the user
            if isempty(val)
                return
            end
            
            %Fully qualified path name or relative?
            [p, ~, ~] = fileparts(val);
            if isempty(p) %Default to current working directory
                val = fullfile(pwd, val);
            end
            
        end
        function val = get.CurrentResults(obj)    %get.CurrentResults
           %get.CurrentResults Retrieves the handle to the current results 
           %are selected 
           
           %Simply index the list in the underlying results viewer object
           val = obj.ResultsView.Results(obj.Selection);
           
        end
        function val = get.AllBeamResults(obj)    %get.AllBeamResults
            %get.AllBeamResults Get method for the dependent property 
            %'AllBeamResults'.
            
            %Grab the beam results from each results set
            val = arrayfun(@(r) r.BeamResults, obj.CurrentResults, 'Unif', false);
            
            %Remove empties
            idx = cellfun(@(x) isempty(x), val);
            val(idx) = [];
            
        end
        function val = get.AllBeams(obj)          %get.AllBeams
            %get.AllBeams Get method for the dependent property 'AllBeams'.
            
            val = cellfun(@(br) [br.Beam], obj.AllBeamResults, 'Unif', false);

        end
        function val = get.AllBeamNames(obj)      %get.AllBeamNames
            %get.AllBeamNames Get method for the dependent property
            %'AllBeamNames'.
            
            val = cellfun(@(x) {x.Name}, obj.AllBeams, 'Unif', false);
        end
        function val = get.CommonBeamNames(obj)   %get.CommonBeamNames
            %get.CommonBeamNames Get method for the dependent property
            %'CommonBeamNames'.
            
            %Start with all the beam names
            allNames = obj.AllBeamNames;
            
            %Find the unique values
            uniqueNames = unique([allNames{:}]);
            
            %Find the common beam names across all beams and return this
            nTypes = cellfun(@(x) numel(intersect(uniqueNames, x)), allNames);
                        
            %Select the cell with the lowest number of interestions
            [~, ind] = min(nTypes);
            val  = allNames{ind};
            
        end
        function val = get.CommonResultTypes(obj) %get.CommonResultTypes
            %get.CommonResultTypes Get method for the dependent property
            %'CommonResultTypes'.
            
            %Grab the name of all BeamProperty types and the unique values
            allTypes = cellfun(@(x) x.BeamPropertyTypes, obj.AllBeamResults, 'Unif', false);
            uniqueTypes = unique(horzcat(allTypes{:}), 'stable');
            
            %Find the common types across all beam results and return this
            nTypes = cellfun(@(x) numel(intersect(uniqueTypes, x)), allTypes);
            
            %Select the cell with the lowest number of interestions
            [~, ind] = min(nTypes);
            val  = allTypes{ind};
            
        end        
        function val = get.ExportAWI(obj)         %get.ExportAWI
            %get.ExportAWI Get method for the dependent property
            %'ExportAWI'.
            
            val = logical(obj.hAWI.Value);
        end
        function val = get.ExportCSV(obj)         %get.ExportCSV
            %get.ExportCSV Get method for the dependent property
            %'ExportCSV'.
            
            val = logical(obj.hCSV.Value);
        end
        function val = get.ExportExcel(obj)       %get.ExportExcel
            %get.ExportExcel Get method for the dependent property
            %'ExportExcel'.
            
            val = logical(obj.hExcel.Value);
        end
        function val = get.ExportFlag(obj)        %get.ExportFlag
            %get.ExportFlag Get method for the dependent property
            %'ExportFlag'.
            
            val = [obj.ExportAWI, obj.ExportCSV, obj.ExportExcel];
            
        end
        function val = get.FileOpenIcon(obj)      %get.IconFileDir
            %get.IconFileDir Get method for the dependent property
            %'IconFileDir'.
            
            %Where are all Matlab icons stored?
            [~, path] = ipticondir;
            
            %Return the image data as an array
            val = imread(fullfile(path, obj.FileOpenIconName));
            
            %Convert to double & normalise by 256^2
            val = double(val) ./ (256^2);
        end
    end
    
    methods % construction
        function obj = ResultsExporter(ResView, varargin)
            
            %Check inputs
            if nargin == 0 || ~isa(ResView, 'awi.view.BeamResultsViewer')
                error(['Cannot initiate an instance of the ', ....
                    '''AWI Results Exporter'' without providing a valid ', ...
                    'instance of ''awi.view.BeamResultsViewer''.']);                
            end
            
            %Grab valid parameter inputs
            
            %Store a reference to the results viewer & honour current
            %selection
            obj.ResultsView = ResView;
            obj.Selection   = obj.ResultsView.Selection;
            
            %Make our own export manager...
            obj.hFig = figure('Name', 'AWI Export Manager', ...
                'WindowStyle', 'normal', ...
                'NumberTitle', 'off'   , ...
                'MenuBar'    , 'none'  , ...
                'ToolBar'    , 'none'  , ...
                'CloseRequestFcn', @obj.cbOnFigClose);
            
            %hBoxFlex for all main components
            hBox = uix.HBoxFlex('Parent', obj.hFig);
            
            % Box for selecting results sets
            vBox = uix.VBox('Parent', hBox, 'Padding', obj.Padding, 'Spacing', obj.Spacing);
            
            %Tell the user what to do!
            uicontrol('Parent', vBox, 'Style', 'text', ...
                'String', 'Select the results sets for exporting', ...
                'HorizontalAlignment', 'left');
            
            %List the results (Allow multiple selection)
            obj.hList = uicontrol('Parent', vBox, 'style', 'listbox', ...
                'String'  , {obj.ResultsView.Results.Name}, ...
                'Max'     , 2, ...
                'Value'   , obj.Selection, ...
                'Callback', @obj.cbOnSelectionChanged);
            
            %Button to quickly select everything
            uicontrol('Parent', vBox, 'Style', 'Push', ...
                'String', 'Select All', ...
                'Callback', @obj.cbSelectAll);
            
            %Buttons for getting/displaying the results directory
            hBox_ = uix.HBox('Parent', vBox, 'Padding', obj.Padding, 'Spacing', obj.Spacing);
            uicontrol('Parent', hBox_, 'Style', 'text', ...
                'String', 'Select a folder for the results', ...
                'HorizontalAlignment', 'left');
            uicontrol('Parent', hBox_, 'Style', 'push',  ...
                'CData'   , obj.FileOpenIcon, ... %Sets the icon as the image
                'Callback', @obj.cbSelectResDir);
            hBox_.Widths = [-1, obj.ButtonHeight];            
            obj.hResDir = uicontrol('Parent', vBox, 'Style', 'edit', ...
                'HorizontalAlignment', 'left', ...
                'BackgroundColor', 'w');            
                        
            %Configure heights
            vBox.Heights = [obj.ButtonHeight, -1, obj.ButtonHeight, obj.ButtonHeight, obj.ButtonHeight];
            
            %Parent box for the directory and export configurables
            vBoxExp = uix.VBox('Parent', hBox, 'Padding', obj.Padding, 'Spacing', obj.Spacing);
            
            %Vertical box for selecting beams & results
            vBoxSel = uix.VBox('Parent', vBoxExp, 'Padding', obj.Padding, 'Spacing', obj.Spacing);
            
            %Box for the beam names (e.g. Fuselage, Wing, HTP, etc.)
            uicontrol('Parent', vBoxSel, 'Style', 'text', ...
                'String', 'Beam Names', ...
                'HorizontalAlignment', 'left');
            obj.hBeamName = uicontrol('Parent', vBoxSel, 'Style', 'list', ...
                'Max'   , 2, ...
                'String', '-Select a results set-');
            
            %Box for the results types (e.g. Internal Loads, Aerodynamics)
            uicontrol('Parent', vBoxSel, 'Style', 'text', ...
                'String', 'Results Types', ...
                'HorizontalAlignment', 'left');
            obj.hResType = uicontrol('Parent', vBoxSel, 'Style', 'list',...
                'Max'   , 2, ...
                'String', '-Select a results set-');
                        
            %Configure heights
            vBoxSel.Heights = [obj.ButtonHeight, -1, obj.ButtonHeight, -1];
            
            %Buttons for configuring the export spec
            vBBox = uix.VButtonBox('Parent', vBoxExp, ...
                'Padding', obj.Padding, ...
                'Spacing', obj.Spacing, ...
                'HorizontalAlignment', 'left');
            obj.hAWI   = uicontrol('Parent', vBBox, 'Style', 'check', 'String', 'AWI Framework (.mat)', 'Enable', 'off');
            obj.hCSV   = uicontrol('Parent', vBBox, 'Style', 'check', 'String', 'CSV (.csv)');
            obj.hExcel = uicontrol('Parent', vBBox, 'Style', 'check', 'String', 'Excel (.xlsx)');
            pos = vBoxExp.Position;
            vBBox.ButtonSize = [pos(3), obj.ButtonHeight];
            
            %Button for passing selection back to the ResultsViewer
            uicontrol('Parent', vBoxExp, 'Style', 'push', ...
                'String', 'Export', ...
                'Callback', @obj.cbConfigureExportables);
            
            %Configure heights
            vBoxExp.Heights = [-1, 4 * obj.ButtonHeight, obj.ButtonHeight];
            
            %Update the GUI
            update(obj);
            
            %Allow the execution to hang until the user has finished
            %interacting with it 
            uiwait(obj.hFig);
        end
    end
    
    methods % update
        function update(obj)
            %update Updates all aspects of the GUI
            
            %Update the list boxes
            obj.hList.Value      = obj.Selection;
            obj.hBeamName.String = obj.CommonBeamNames;
            obj.hResType.String  = obj.CommonResultTypes;
            
        end
    end
    
    methods % callbacks
       
        function cbOnSelectionChanged(obj, ~, ~)
            %cbOnSelectionChanged Triggers the update of the other list
            %boxes so that the correct data is displayed.
            
            %Update the selection 
            obj.Selection = obj.hList.Value;
            
            %Force an update
            update(obj);
        end
        
        function cbSelectAll(obj, ~, ~) 
            %cbSelectAll Callback for the 'Select All' button 
            
            %Simple!            
            obj.hList.Value = 1 : numel(obj.hList.String);
            obj.Selection   = obj.hList.Value;
            
            %Force an update
            update(obj);
        end
        
        function cbSelectResDir(obj, ~, ~)
            %cbSelectResDir Opens up a dialogue box for selecting the
            %directory where the results will be exported to.
            
            %Request a directory for the results to be sent to
            path = uigetdir(pwd, ['Select a folder for the results ', ...
                'to be exported to.']);
            
            %Escape route
            if path == 0
                return
            end
            
            %Propogate the change to the uicontrol object that displays the
            %path
            obj.hResDir.String = path;
            
        end
        
        function cbConfigureExportables(obj, ~, ~)
            %cbConfigureExportables Passes the selection from the
            %'ResultsExporter' to the 'ResultsViewer' for carrying out the
            %export and closes the 'ResultsExporter' to trigger 'uiwait'.
            
            %Pass information back the the 'ResultsViewer'
            if isempty(obj.ResDir)           %Results directory provided?
                errordlg(['You must select a directory for the results ', ...
                    'to be exported to before you can export the data.'], ...
                    'No results directory selected');
                return
            end
            if exist(obj.ResDir, 'dir') ~= 7 %Results directory exist?
                %Ask the uicontrol object is of style 'edit' we need to 
                %check the user hasn't supplied a folder that doesn't 
                %actually exist!
                errordlg(sprintf(['The folder ''%s'' does not exist. ' , ...
                    'Please provide a valid folder for the results to ', ...
                    'be exported to'], obj.ResDir),  ...
                    'Results folder does not exist');
                return
            end
            if ~any(obj.ExportFlag)          %Export spec defined?
                errordlg(['You must select a file type for the results ', ...
                    'to be exported to before you can export the data.'], ...
                    'No file type selected');
                return
            end
            
            %Which results sets?
            obj.ResultsView.Res2Export = obj.CurrentResults;
            
            %Which beams?
            obj.ResultsView.BeamNames2Export = obj.hBeamName.String(obj.hBeamName.Value);
            
            %Which types of results?
            obj.ResultsView.ResTypes2Export = obj.hResType.String(obj.hResType.Value);
            
            %What file format?
            obj.ResultsView.ExportAWI   = obj.ExportAWI;
            obj.ResultsView.ExportCSV   = obj.ExportCSV;    
            obj.ResultsView.ExportExcel = obj.ExportExcel;

            %Where should the data be saved?
            obj.ResultsView.ExportPath = obj.ResDir;
            
            %Delete the Export Manager and return to the Results Viewer
            %   - Do not want to call 'close' as this will trigger the
            %   'CloseRequestFcn'.
            delete(obj.hFig);   
            
        end
        
        function cbOnFigClose(obj, ~, ~)
            %cbOnFigClose Ensures that no export is carried out if the user
            %closes the Export Manager.
            
            %Make sure nothing gets exported!
            obj.ResultsView.Res2Export       = [];
            obj.ResultsView.BeamNames2Export = {};
            obj.ResultsView.ResTypes2Export  = {};
            obj.ResultsView.ExportAWI        = false;
            obj.ResultsView.ExportCSV        = false;    
            obj.ResultsView.ExportExcel      = false;
            obj.ResultsView.ExportPath       = '';
            
            %Close the figure
            delete(obj.hFig);
            
        end
    end
    
end

