classdef SelfTest < matlab.unittest.TestCase & matlab.mixin.SetGet
    %
    % Implements a suite of tests of the mvc class library
    
    properties
        ObjectUnderTest;
        SaveFolder;
        DataFolder;
        ViewManager = @mvc.view.ViewManager;
        
    end
    
    properties (TestParameter, SetAccess = protected)
        
        %Importable files
        ImportFile = {'x.xml'};
                    
    end
        
    methods % construction
        
        function obj = SelfTest(model, varargin)
            
            %Caller may pass in the thing to be tested
            if nargin < 1 || isempty(model)
                
                %Create it
                model = mvc.model.Application;
                
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
                obj.DataFolder = fullfile(fileparts(fileparts(fileparts(fileparts(mfilename('fullpath'))))), 'Data');
            end
            assert(isdir(obj.DataFolder), ['data folder ''', obj.DataFolder, ''' not found']);
            
            %If caller does not want it back
            if nargout == 0
              
                %Show us save folder content (so can see saved files as they appear)
                winopen(obj.SaveFolder);
                
                %What else would you want to do ?
                testResults = run(obj); %#ok<NASGU>
                
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
            if isa(obj.ObjectUnderTest, 'mvc.mixin.Collectable') && obj.ObjectUnderTest.HasChildren
                
                %Yes - start a new session
                obj.ObjectUnderTest = new(obj.ObjectUnderTest);
            
            end
            
        end
        
    end
    
    methods (Test)
        
        function vw = viewThenAddOneOfEverything(obj)
            
            %Good filename ?
            fn = fullfile(obj.SaveFolder, 'viewThenAddOneOfEverything');
            
            %Start from clean sheet
            obj.cleanSession;
            
            %Create the view
            vw = obj.ViewManager(obj.ObjectUnderTest);
            
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
            vw = obj.ViewManager(obj.ObjectUnderTest);
            
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
            vw = obj.ViewManager(obj.ObjectUnderTest);
            
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
            vw = obj.ViewManager(obj.ObjectUnderTest);
            
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
        
%         function vw = viewThenImportThenAnalyse(obj, AnalyseFile)
%             
%             %Good filename ?
%             fn = fullfile(obj.SaveFolder, ['viewThenImportThenAnalyse-', AnalyseFile]);
%             
%             %Start from clean sheet
%             obj.cleanSession;
%             
%             %Create the view
%             vw = obj.ViewManager(obj.ObjectUnderTest);
%             
%             %Import the specified file
%             obj.ObjectUnderTest.import(fullfile(obj.DataFolder, AnalyseFile));
%             
%             %Pass it on
%             obj.doAnalyses(obj.ObjectUnderTest, vw);
%         
%             %Save session to file
%             obj.saveCleanAndReload(obj.ObjectUnderTest, fn);
%          
%             %Export session to file
%             obj.exportCleanAndReimport(obj.ObjectUnderTest, fn);
%             
%             %Exercise the view(s)
%             obj.addOneOfEveryView(vw);
%                                  
%             %If caller does not want the view back
%             if nargout == 0
%                 
%                 %Close it, with force flag to suppress confirmation prompt
%                 vw.close(true);
%                 
%             end
%                         
%         end
        
    end
    
    methods (Static)
        
        function addOneOfEveryView(obj)
            
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
            
            %But don't do everything - just choose a selection at random
            for i = randi(numel(M), 1, 5)
                obj.Selection = M(i);
                drawnow
            end
            
        end
        
        function addOneOfEverything(obj, varargin)
            
            %What is addable ?
            spec = obj.CollectionSpec;
            
            %Only interested in those that are not hidden (I think ?)
            spec([spec.Hidden]) = [];
            
            %To avoid infinite recursion, ignore anything we've already added
            spec(ismember({spec.Description}, varargin)) = [];
            
            %For each in turn
            for i = 1:numel(spec)
            
                %Add one
                item = obj.add(spec(i).Description);
                
                %If this item is itself NOT a leaf node
                if ~item.IsLeafNode
                    
                    %Recurse
                    mvc.model.SelfTest.addOneOfEverything(item, varargin{:}, spec(i).Description);
        
                end
                
            end
                        
        end
        
        function fn = saveCleanAndReload(obj, fn)
            
            %Strip extension from file
            [p,n] = fileparts(fn);
            
            %Save session
            save(obj, fullfile(p,n));
          
            %Take a snapshot
            str = mvc.model.SelfTest.snapshot(obj);
            
            %Recover the actually-used filename
            fn = obj.SerializeFullFile;
            
            %Clean
            obj.Children = [];
            
            %Reload session, passing in the "force" option to suppress prompt
            obj.open(fn, true);
          
            %Compare
            mvc.model.SelfTest.snapshot(obj, str);
            
        end
        
        function fx = exportCleanAndReimport(obj, fn)
            
            %Take a snapshot before we do anything
            str = mvc.model.SelfTest.snapshot(obj);
            
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
                mvc.model.SelfTest.snapshot(obj, str);
            
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
                
                %Limit ourselves to char (for time being)
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
