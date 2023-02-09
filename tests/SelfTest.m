classdef SelfTest < matlab.unittest.TestCase & matlab.mixin.SetGet
    %
    % Implements a suite of tests of the AWI Framework
    % Original code by Phil Rottier at MATHWORKS - Need to update tests
    % significantly.
    %
    % TODO:
    %   - Add a test where we initiate each 'model' class and check it has
    %     a default build.
    %   - Attempt to run each available public method for the 'model'
    %    classes.
    
    properties
        ObjectUnderTest;
        SaveFolder;
        DataFolder;
    end
    
    properties (TestParameter)
        
        %Standard test aircraft defined in XML files
        ImportFile = { ... 'harten_v2.xml', ...     
            ['\\ads\filestore\Engineering\Research\Projects\AWI\', ...
            'Analysis_Tools\UoB_Framework\Models\HARTEN_v2_VLM\' , ...
            '06_weights\FAME\v3_newCG_LC.fame-w.4.00.m1-1-pv8\'  , ...
            '__mirror\newCG_LC.fm4'], ...
            'BUG.xml'      , ...
            'sugar_volt_765-095_RevD_v4.xml'};
        
        %Files that can be analysed
        AnalyseFile = { ...
            ['\\ads\filestore\Engineering\Research\Projects\AWI\', ...
            'Analysis_Tools\UoB_Framework\Models\HARTEN_v2_VLM\' , ...
            '06_weights\FAME\v3_newCG_LC.fame-w.4.00.m1-1-pv8\'  , ...
            '__mirror\newCG_LC.fm4']};
            
    end
        
    methods % construction
        
        function obj = SelfTest(model, varargin)
            
            %Caller may pass in the thing to be tested
            if nargin < 1 || isempty(model)
                
                %Create it
                model = awi.model.Framework;
                
            end
            
            %Assignments ?
            if ~isempty(varargin)
                set(obj, varargin{:});
            end
            
            %Store it internally
            obj.ObjectUnderTest = model;
            
            %Create a folder into which results are saved
            tok = strsplit(class(obj), '.');
            if isempty(obj.SaveFolder)
                obj.SaveFolder = fullfile(pwd, tok{end}, datestr(now, 30));
            end
            mkdir(obj.SaveFolder);
            
            %Where does the data live ?
            if isempty(obj.DataFolder)
                obj.DataFolder = '\\ads\filestore\Engineering\Research\Projects\AWI\Analysis_Tools\UoB_Framework\Models';
                %obj.DataFolder = fullfile(fileparts(fileparts(fileparts(fileparts(mfilename('fullpath'))))), 'Data');
            end
            assert(isdir(obj.DataFolder), ['data folder ''', obj.DataFolder, ''' not found']);
            
            %If caller does not want it back
            if nargout == 0
              
                %Show us save folder content (so can see saved files as they appear)
                winopen(obj.SaveFolder);
                
                %What else would you want to do but run at least one test ?
                mc = metaclass(obj);
                tdx = find(arrayfun(@(x)isprop(x, 'Test') && x.Test, mc.MethodList));
                tests = {mc.MethodList(tdx).Name};
                sel = obj.ObjectUnderTest.listdlg(tests, 'Choose test to run, or CANCEL to run ALL', ...
                    'SelectionMode', 'single');
                
                %Profile ?
                if obj.ObjectUnderTest.confirm('Profile ?')
                    fcn = @profile;
                else
                    fcn = @run;
                end
                
                %Run everything ?
                if isempty(sel)
                    
                    %Yes
                    testResults = fcn(obj); %#ok<NASGU>
                    
                else
                    
                    %Just the one
                    testResults = fcn(obj, mc.MethodList(tdx(sel)).Name); %#ok<NASGU>
                    
                end
                
                %Save test results in results folder
                save(fullfile(obj.SaveFolder, 'testResults'), 'testResults');
                
            end
            
        end
        
        function varargout = profile(obj, varargin)
            
            %Turn on profiler
            profile('on');
            
            %Run whatever
            [varargout{1:nargout}] = run(obj, varargin{:});
            
            %What happened ?
            profile('off');
            profile('report');
            
            %And save the profile info to the results folder
            profileInfo = profile('info'); %#ok<NASGU>
            save(fullfile(obj.SaveFolder, 'profileInfo'), 'profileInfo');
            
        end
        
    end
        
    methods % (Access = protected) % helpers
        
        function cleanSession(obj)

            %Anything to worry about ?
            if obj.ObjectUnderTest.HasChildren
                
                %Yes - start a new session
                obj.ObjectUnderTest = new(obj.ObjectUnderTest);
            
            end
            
        end
        
    end
    
    methods (Test)
        
        function vw = viewThenAddOneOfEverything(obj)
            %viewThenAddOneOfEverything Unit test No. 1
            %
            % The following operations are performed as part of this test:
            %   - Initialise an AWI session
            %   - Add one of every type of view that is supported by AWI
            %   - Add one of every object that is supported by AWI
            %   - Save the session to the .awi format and attempt to reload
            %   - Export the session to all available formats and import
            %     those files back into the framework.
            %
            % N.B. This unit test will take longer than its 'sister' test
            % 'addOneOfEverythingThenView' as the GUI elements are present
            % whilst objects are added to the collection. This has a
            % significant overhead as the event-listener methods in the MVC
            % package are poorly optimised.
            
            %Good filename ?
            fn = fullfile(obj.SaveFolder, 'viewThenAddOneOfEverything');
            
            %Start from clean sheet
            obj.cleanSession;
            
            %Create the view
            vw = awi.view.Framework(obj.ObjectUnderTest);
            
            %Exercise the view(s), before we add any content
            obj.addOneOfEveryView(vw);
            
            %Pass it on
            obj.addOneOfEverything(obj.ObjectUnderTest);
        
            %Save session to file
            obj.saveCleanAndReload(obj.ObjectUnderTest, fn);
         
            %Export session to file
            obj.exportCleanAndReimport(obj.ObjectUnderTest, fn);
                     
            %If caller does not want the view back
            if nargout == 0
                
                %Close it, with force flag to suppress confirmation prompt
                vw.close(true);
                
            end
            
        end
        
        function vw = addOneOfEverythingThenView(obj)
            %viewThenAddOneOfEverything Unit test No. 2
            %
            % The following operations are performed as part of this test:
            %   - Initialise an AWI session
            %   - Add one of every object that is supported by AWI     
            %   - Save the session to the .awi format and attempt to reload
            %   - Export the session to all available formats and import
            %     those files back into the framework.
            %   - Add one of every type of view that is supported by AWI
            
            %Good filename ?
            fn = fullfile(obj.SaveFolder, 'addOneOfEverythingThenView');
            
            %Start from clean sheet
            obj.cleanSession;
            
            %Pass it on
            obj.addOneOfEverything(obj.ObjectUnderTest);
         
            %Save session to file
            obj.saveCleanAndReload(obj.ObjectUnderTest, fn);
         
            %Export session to file
            obj.exportCleanAndReimport(obj.ObjectUnderTest, fn);
                     
            %Create the view
            vw = awi.view.Framework(obj.ObjectUnderTest);
            
            %Exercise the view(s)
            obj.addOneOfEveryView(vw);
                        
            %If caller does not want the view back
            if nargout == 0
                
                %Need to pause else you won't see it
                pause(1);
                
                %Close it, with force flag to suppress confirmation prompt
                vw.close(true);
                
            end
          
        end
        
        function vw = viewThenImport(obj, ImportFile)
            
            %Good filename ?
            fn = fullfile(obj.SaveFolder, ['viewThenImport-', ImportFile]);
            
            %Start from clean sheet
            obj.cleanSession;
            
            %Create the view
            vw = awi.view.Framework(obj.ObjectUnderTest);
            
            %Import the specified file
            obj.ObjectUnderTest.import(fullfile(obj.DataFolder, ImportFile));
        
            %Save session to file
            obj.saveCleanAndReload(obj.ObjectUnderTest, fn);
         
            %Export session to file
            obj.exportCleanAndReimport(obj.ObjectUnderTest, fn);
            
            %Exercise the view(s)
            obj.addOneOfEveryView(vw);
                                 
            %If caller does not want the view back
            if nargout == 0
                
                %Close it, with force flag to suppress confirmation prompt
                vw.close(true);
                
            end
                        
        end
        
        function vw = importThenView(obj, ImportFile)
            
            %Good filename ?
            fn = fullfile(obj.SaveFolder, ['importThenView-', ImportFile]);
            
            %Start from clean sheet
            obj.cleanSession;
            
            %Import the specified file
            obj.ObjectUnderTest.import(fullfile(obj.DataFolder, ImportFile));
        
            %Save session to file
            obj.saveCleanAndReload(obj.ObjectUnderTest, fn);
         
            %Export session to file
            obj.exportCleanAndReimport(obj.ObjectUnderTest, fn);
            
            %Create the view
            vw = awi.view.Framework(obj.ObjectUnderTest);
            
            %Exercise the view(s)
            obj.addOneOfEveryView(vw);
                                 
            %If caller does not want the view back
            if nargout == 0
                
                %Need to pause else you won't see it
                pause(1);
                
                %Close it, with force flag to suppress confirmation prompt
                vw.close(true);
                
            end
          
        end
        
        function vw = viewThenImportThenAnalyse(obj, AnalyseFile)
            
            %Good filename ?
            fn = fullfile(obj.SaveFolder, ['viewThenImportThenAnalyse-', AnalyseFile]);
            
            %Start from clean sheet
            obj.cleanSession;
            
            %Create the view
            vw = awi.view.Framework(obj.ObjectUnderTest);
            
            %Import the specified file
            if isempty(fileparts(AnalyseFile))
                obj.ObjectUnderTest.import(fullfile(obj.DataFolder, AnalyseFile)); %requires file to be in Data folder
            else
                obj.ObjectUnderTest.import(AnalyseFile); % allow "which" to locate file
            end
            
            %Pass it on
            obj.doAnalyses(obj.ObjectUnderTest, vw, fn);
        
            %Save session to file
            obj.saveCleanAndReload(obj.ObjectUnderTest, fn);
         
            %Export session to file
            obj.exportCleanAndReimport(obj.ObjectUnderTest, fn);
            
            %Exercise the view(s)
            obj.addOneOfEveryView(vw);
                                 
            %If caller does not want the view back
            if nargout == 0
                
                %Close it, with force flag to suppress confirmation prompt
                vw.close(true);
                
            end
                        
        end
        
        function vw = viewThenImportThenViewResults(obj, ImportFile)
            
            %Good filename ?
            fn = fullfile(obj.SaveFolder, ['viewThenImportThenViewResults-', ImportFile]);
            
            %Start from clean sheet
            obj.cleanSession;
            
            %Create the view
            vw = awi.view.Framework(obj.ObjectUnderTest);
            
            %Import the specified file
            if isempty(fileparts(ImportFile))
                obj.ObjectUnderTest.import(fullfile(obj.DataFolder, ImportFile)); %requires file to be in Data folder
            else
                obj.ObjectUnderTest.import(ImportFile); % allow "which" to locate file
            end
            
            %Pass it on
            %obj.doAnalyses(obj.ObjectUnderTest, vw, fn);
        
            %Save session to file
            %obj.saveCleanAndReload(obj.ObjectUnderTest, fn);
         
            %Export session to file
            %obj.exportCleanAndReimport(obj.ObjectUnderTest, fn);
            
            %Exercise the view(s)
            %obj.addOneOfEveryView(vw);
                                 
            %If caller does not want the view back
            if nargout == 0
                
                %Close it, with force flag to suppress confirmation prompt
                vw.close(true);
                
            end
                        
        end
        
    end
    
    methods (Static)
        
        function viewResults(obj, vw, fn)
            %viewResults Static method for viewing all of the available
            %results in the AWI session and exercising the full
            %functionality of the results views.
            %
            % The following views are available to viewing results:
            %   - 'awi.view.BeamResultsViewer'
            
                        
        end
        
        function doAnalyses(obj, vw, fn, varargin)
           
            %Add a mass analysis
            ma = obj.massAnalysis(varargin{:});

            %Do it - actually there's nothing to do
            % ma.analyse;
            
            %Creating the Size view
            sz = obj.sizeAnalysis(varargin{:});
            
            %Do it
            sz.analyse;
            
            %Create the Trim view
            tm = obj.trimAnalysis(varargin{:});
            
            %Do it
            tm.analyse;
            
            %Add a results viewer
            rv = vw.addView('Results Viewer');
            
            %Select all results
            rv.Selection = rv.Results;
            
            %Name of Word doc to copy view to
            doc = strrep(fn, '.', '_');
            
            %Create one of each arrangement
            for i = 1:numel(rv.SupportedArrangements)
                
                %Assign the arrangement
                rv.Arrangement = i;
                
                %Copy view to Word
                doc = vw.copyView(rv, 'Word', doc, ['Arrangement ''', rv.Arrangement.Name, '''']);
                
                %Export view to Excel
                b(i) = rv.export(strrep(fn, '.', '_'), 'Sheet', rv.Arrangement.Name);
                
            end
            
            %Close document
            vw.copyView([], 'Word', doc);
            
            %Ensure success
            assert(all(b), 'Excel export failed');
            
        end
        
        function addOneOfEveryView(obj)
            %addOneOfEveryView Static method for adding a single instance
            %of every supported view to the AWI session.
            %   
            % Steps performed: 
            %   - Grab supported views.
            %   - Add one of every view to the session whilst exercising
            %     the ability to move views around the arrangement.
            %   - Test the drawing of elements of the Framework.
            
            %What views are available ?
            vws = obj.SupportedViews;
            
            %Options for position of views ?
            pos = {'right', 'below', ''};
            ipos = 1;
            
            %Add one of each
            for i = 1:size(vws,1)
                
                %Add the view
                vw = obj.addView(vws{i,1});
                
                %Move it somewhere
                if ~isempty(pos{ipos})
                    obj.moveView(vw, pos{ipos});
                end
                
                %Clock-increment the position counter
                ipos = ipos + 1;
                if ipos > numel(pos)
                    ipos = 1;
                end
                
            end
            
            %Verify that we can programatically select members of the model
            M = obj.Model.flatlist;
            
            %Down-select to items that actually appear in tree
            M(~[M.DrawChildrenInTree]) = [];
            
            %But don't do everything - just choose a selection at random
            for i = randi(numel(M), 1, 5)
                obj.Selection = M(i);
                drawnow
            end
            
        end
        
        function created = createOneOfEverything(obj, varargin)
            
            %We've already created one of these
            created{1} = class(obj);
            
            %Then a placeholder for a list of additional things created
            created{2} = {};
            
            %What is addable ?
            spec = obj.CollectionSpec;
            
            %To avoid infinite recursion, ignore anything we've already created
            
            
        end
        
        function added = addOneOfEverything(obj, varargin)
            %addOneOfEverything Static method for adding one of every
            %collectable object to the AWI session.
            %
            % Steps performed:
            %   - Grab collection specification of the current object.
            %   - Check if we have already added one of the items in the
            %     specification and if yes, remove so we don't add again.
            %   - For each remaing item in the specification, add it to the
            %     session and call 'get' to check that we can successfully
            %     retrieve a valid object.  
            %   - If we can add more than one of the items then make a
            %     duplicate of that item and remove the original. 
            %   - Duplicate the item again and call 'get'. Compare this
            %     with the original call to 'get' to ensure we have made a
            %     deep copy of the object.
            %   - Recursively call this function for every item in the
            %     collection spec.            
            
            %Return arg starts with name of this class
            added{1} = class(obj);
            
            %Then a placeholder for a list of things actually added
            added{2} = {};
            
            %What is addable ?
            spec = obj.CollectionSpec;
            
            %Only interested in those that are not hidden (I think ?)
            spec([spec.Hidden]) = [];
            
            %To avoid infinite recursion, ignore anything we've already added
            % to an object of this class already
            for i = find(strcmp(class(obj), varargin(1:2:end)))           
                spec(ismember({spec.Description}, varargin{i * 2})) = [];
            end
            
            %For each in turn of whatever is left
            for i = 1:numel(spec)
            
                %What's happening ?
                disp([added{1}, ' adding ', spec(i).Description, '...']);
                
                %Add one
                item = obj.add(spec(i).Description);
                
                %Ensure we can call 'get' on this thing
                X = get(item);
                
                %Make a note of what was added
                added{2}{end+1} = spec(i).Description;
                
                %Can we add more than one of these things ?
                if spec(i).MaxNumber > 1
                    
                    %Verify that we can duplicate it
                    dupitem = item.duplicateThis;
                    
                    %Verify that we can delete the original
                    item.removeThis(true);
                                        
                    %And replace it with the duplicate
                    % (therefore ensuring that we really have made a deep copy)
                    item = dupitem;
                    
                    %Call get again
                    Y = get(item);
                    
                    %Verify that result of original and subsequent call to 'get' are same
                    if ~isequaln(X, Y)
                                                
                        %For objects that aggregate other objects we would
                        %not neccessarily expect 'isequaln' to return true
                        %as the duplicate and the original generate their
                        %own (unique) objects upon construction (e.g. with
                        %the ControlSurface class - and it's subclasses!)
                        %
                        %   - Check that this is the case
                        
                        %Test each field of the structures 'X' & 'Y'
                        fnames = fieldnames(X);
                        tf     = false(size(fnames));
                        for iF = 1 : numel(fnames)
                            tf(iF) = isequal(X.(fnames{iF}), Y.(fnames{iF}));
                        end
                        
                        %What is not equal?
                        fnames = fnames(~tf);
                        val    = cell(numel(fnames), 2);
                        tf     = false(size(fnames));
                        for iF = 1 : numel(fnames)
                            %Index the structures
                            val{iF, 1} = X.(fnames{iF});
                            val{iF, 2} = Y.(fnames{iF});
                            %If they are MATLAB objects then disregard,
                            %otherwise flag it up.
                            tf(iF) = or(~isobject(val{iF, 1}), ~isobject(val{iF, 2})); %TODO - This is still not good enough. Ideally we would have a list of allowable objects and properties that can skip this assertion.
                        end
                        
                        if any(tf)
                            %Ideally would assert, but need to work on it a bit more
                            %   - TODO : Use tf to list all incorrect
                            %   properties
                            warning([class(item), ': duplication error']);
                        end
                        
                    end
                    
                end
                
                %If this item is itself NOT a leaf node
                if ~item.IsLeafNode
                    
                    %Recurse
                    also_added = awi.model.SelfTest.addOneOfEverything(item, varargin{:}, added{:});
        
                    %Make a note of what happened (if anything)
                    if ~isempty(also_added{2})
                        added = [added, also_added];
                    end
                    
                end
                
            end
                        
        end
        
        function fn = saveCleanAndReload(obj, fn)
            %saveCleanAndReload Static method for saving the AWI session to
            %a '.mat' file, removing the children of the current session
            %and then reimporting the session data for comparison.
            %
            % Actions performed:
            %   - Save the session to a '.mat' file
            %   - Take a snapshot of the session
            %   - Remove all children of the session
            %   - Reimport the session data from the '.mat' file
            %   - Take a snapshot and compare with the previous instance.
            %
            % see also: mvc.mixin.Serializable, awi.model.SelfTest.snapshot
            %           
            
            %Strip extension from file
            [p,n] = fileparts(fn);
            
            %Save session
            save(obj, fullfile(p,n));
          
            %Take a snapshot
            str = awi.model.SelfTest.snapshot(obj);
            
            %Recover the actually-used filename
            fn = obj.SerializeFullFile;
            
            %Clean
            obj.Children = [];
            
            %Reload session, passing in the "force" option to suppress prompt
            obj.open(fn, true); %TODO - This line causes several warnings about heterogeneous arrays and prop/event listeners. Need to fix.
          
            %Compare
            awi.model.SelfTest.snapshot(obj, str);
            
        end
        
        function fx = exportCleanAndReimport(obj, fn)
            
            %Take a snapshot before we do anything
            str = awi.model.SelfTest.snapshot(obj);
            
            %Strip extension from file
            [p,n] = fileparts(fn);
            
            %What export file extensions are supported ?
            ex = obj.ExportFileExtensions;
            
            %Placeholder
            fx = cell(size(ex));
            
            %For each in turn
            for i = 1:numel(ex)
                
                %Export session
                export(obj, fullfile(p, [n, ex{i}]));
                
                %Recover the actually-used filename
                fx{i} = obj.ExportFullFile;
                
            end
            
            %Make a note of the original content
            C = detach(obj.Children);
            
            %Which import file extensions are supported ?
            ei = obj.ImportFileExtensions;
            
            %For file formats that are both exportable and importable
            [e, edx] = intersect(ex, ei, 'stable');
            for i = 1:numel(e)
                
                %Clean sheet
                obj.Children = [];
                
                %Reimport session, passing in the "force" option to suppress prompt
                args = {'force', true};
                
                %And the 'all' option for sheets, if we are importing Excel
                if strncmpi(e{i}, '.xls', 4)
                    args(end+1:end+2) = {'Sheet', 'all'};
                end
                
                %Pass it on
                obj.import(fx{edx(i)}, args{:});
                
                %Compare
                awi.model.SelfTest.snapshot(obj, str);
            
            end
            
            %Putback the originals
            obj.Children = C;
            
        end
        
        function [str, diff] = snapshot(obj, ref, prefix)
            %
            %Return a cellstr of textual descriptions of items in collection.
            % Useful for comparing content before and after save/reload
            % or export/reimport, whilst avoiding the need to retain a copy
            % of the actual object
            
            %Get collection as flat list
            lst = obj.flatlist(2);
            
            %Take a snapshot of each
            str = arrayfun(@i_snapshot, lst, 'UniformOutput', false);
            
            %Not bothered about any empties
            str(cellfun(@isempty, str)) = [];
            
            %If a reference snapshot has been provided
            if nargin > 1
                
                %Compare
                if isequaln(str, ref)
                    
                    %Nothing to report
                    diff = {};
                    
                else
                    
                    %Build a diff report
                    b = ~ismember(str, ref);
                    diff = {['Number of expected items: ', num2str(numel(ref))]};
                    diff{end+1} = ['Number of actual items: ', num2str(numel(str))];
                    diff{end+1} = ['Number of actuals present in expected: ', num2str(sum(~b))];
                    diff{end+1} = ['Number of actuals not present in expected: ', num2str(sum(b))];
                    if any(b)
                        diff(end+1:end+sum(b)) = cellfun(@(x)[' ', x], str(b), 'UniformOutput', false);
                    end
                    b = ~ismember(ref, str);
                    diff{end+1} = ['Number of expected not present in actual: ', num2str(sum(b))];
                    if any(b)
                        diff(end+1:end+sum(b)) = cellfun(@(x)[' ', x], ref(b), 'UniformOutput', false);
                    end
                    
                    %If caller did not want difference details back
                    if nargout < 2
                        
                        %Create a message
                        msg = sprintf('%s\n', diff{:});
                        
                        %Caller may provide a prefix
                        if nargin > 2
                            msg = [prefix, ': ', msg];
                        end
                        
                        %Turn this into an error in due course
                        warning(msg);
                    end
                    
                end
                
            end
            
            function str = i_snapshot(obj)
            
                %What properties to snap ?
                prp = allProperties(obj);
                val = get(obj, prp);
                
                %Find all logicals and convert to numerics
                b = cellfun(@islogical, val);                
                val(b) = cellfun(@(x) double(x), val(b), 'Unif', false);
                
                %Find all numerics and convert to char
                b = cellfun(@isnumeric, val);
                val(b) = cellfun(@(x) num2str(x), val(b), 'Unif', false);
                
                %Limit ourselves to char (for time being) - TODO Add all
                %other data types (table, structure, objects, cells, etc.)
                b = cellfun(@ischar, val);
                
                %Down-select
                prp = prp(b);
                val = val(b);
                
                %Hence a snapshot string would be
                pav = [prp; val];
                str = sprintf('%s = %s; ', pav{:});
                
            end
            
        end
        
    end
    
end
