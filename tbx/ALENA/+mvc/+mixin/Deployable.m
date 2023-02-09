classdef (ConstructOnLoad) Deployable < matlab.mixin.SetGet
    %
    % Deployable provides functionality useful for application deployment:
    %
    
    properties (Transient)
    
        %Need a simple entry point function
        EntryPointFcn;
        
        %Deployment targets - reduce this list in application-specific subclasses, if required
        DeploymentTargets = {'MATLAB App', 'MATLAB Toolbox', 'Compiler Project', 'Executable', ...
            'm code', 'p code'};
        
    end
    
    methods % construction / destruction
        
        function obj = Deployable(varargin)
            
            %Extend context (if applicable)
            if isa(obj, 'mvc.mixin.Contextable')
                addContext(obj, Inf, ... % Specify a super-low priority, so this appears at end
                    '|Deploy>Settings...', 'deploymentSettings', ...
                    '|Deploy>Deploy...', 'deploy');
            end
            
            %If serializable
            if isa(obj, 'mvc.mixin.Serializable')
                
                %Ensure EntryPointFcn is propagated into any new instances created
                obj.CopyOnNew{end+1} = 'EntryPointFcn';
                
            end
            
        end
        
    end
    
    methods (Sealed)
        
        function deploymentSettings(obj)
            
            %There is more to it than this, but for time being
            str = {'Entry point function name:'};
            prp = {'EntryPointFcn'};
            val = get(obj, prp);
            newval = inputdlg(obj, val, str);
            if isempty(newval)
                return;
            end
            
            %What has changed ?
            bSame = cellfun(@isequal, val, newval);
            if all(bSame)
                return;
            end
            
            %Apply changes
            pav = [prp; newval];
            pav(:,bSame) = [];
            set(obj, pav{:});
            
        end
        
        function deploy(obj, target, dd, varargin)
            
            %Need a simple entry point function
            assert(~isempty(obj.EntryPointFcn), 'entry point function not set');
            assert(exist(obj.EntryPointFcn, 'file') == 2, 'entry point function not found');
            
            %Deploy as what ?
            if nargin < 2 || isempty(target)
                
                %Ask the user
                target = choose(obj, 'Deployment target...', obj.DeploymentTargets{:});
                
                %Cancelled ?
                if isempty(target)
                    return;
                end
            
            else
                
                %Look for match
                idx = find(~cellfun(@isempty, strfind(lower(obj.DeploymentTargets), lower(target))));
                
                %Valid
                assert(~isempty(idx), 'specified target not recognised');
                assert(numel(idx) == 1, 'specified target not unambiguously recognised');
                
                %So the target is
                target = obj.DeploymentTargets{idx};
                
            end

            %Ensure target is specified in the form of a deployment function
            if ischar(target)
                
                %Map to corresponding deployment function
                fcn = str2func(matlab.lang.makeValidName(['deploy to ', target]));
                
            elseif isa(target, 'function_handle')
                
                %Nothing to do
                fcn = target;
                
            else
                %Nothing else will do
                error('bad input - target must be either a char or a function handle');
            end
            
            %Recycle deployment folder
            persistent DD;
            if isempty(DD)
                DD = pwd;
            end
            
            %Caller specifid deployment directory ?
            if nargin < 3 || isempty(dd)
                
                %No - ask the user
                dd = uigetdir(DD, 'Deployment directory');
                
                %Cancelled ?
                if isempty(dd) || isequal(dd, 0)
                    return;
                end
            
            else
                assert(isdir(dd), 'deployment directory not found');
            end
            
            %Make a note for next time
            DD = dd;
            
            %Create a unique deployment directory
            dd = fullfile(dd, datestr(now, 30));
            if ~mkdir(dd)
                error(['unable to create deployment directory ''', dd, '''']);
            end
            
            %Do something with it
            feval(fcn, obj, dd, varargin{:});
            
            %Show us the result
            winopen(dd);
            
        end

        function deployToMCode(obj, dd, varargin)
        
            %Where is the entry point function located ?
            src = fileparts(which(obj.EntryPointFcn));
            
            %Start from a complete copy of the code folder in the deployment location
            [bOK, msg] = copyfile(src, dd);
            assert(bOK, ['Failed to copy source to destination folder: ', msg]);
        
            %But this is likely to include many classes (e.g. from mvc.mixin .view packages)
            % that are not actually required by whatever app is being deployed,
            % so get a list of all files in the depoyment folder
            fdd = i_dir(dd);            
            function fn = i_dir(src, fn)
                
                %Anything yet ?
                if nargin < 2
                    
                    %No
                    fn = {};
                    
                end
                
                %What have we got ?
                f = dir(src);
                
                %Unpick
                for i = 1:numel(f)
                    
                    %If an actual file
                    if ~f(i).isdir
                        
                        %Ignore asv files
                        [~,~,e] = fileparts(f(i).name);
                        if strcmpi(e, '.asv')
                            continue;
                        end
                        
                        %Fully-qualified name
                        fn{end+1} = fullfile(src, f(i).name);
                        
                    elseif ~ismember(f(i).name, {'.', '..'})
                        
                        %Recurse
                        fn = i_dir(fullfile(src, f(i).name), fn);
                        
                    end
                    
                end
                
            end
            
            %Not interested in anythin in +test package
            testFolders = unique(cellfun(@fileparts, fdd(strncmp(fullfile(dd, '+test'), fdd, numel(fullfile(dd, '+test')))), 'UniformOutput', false));
            cellfun(@(x)rmdir(x, 's'), testFolders);
            
            %Get list of classes associated with application
            cls = superclassList(obj);
            
            %Not interested in Debuggable or Deployable
            cls(ismember(cls, {'mvc.mixin.Debugable', 'mvc.mixin.Deployable'})) = [];
            
            %TODO - get addables (applicable to Collectables only)
            
            %Find a view manager associated with this object
            mgr = viewManager(obj);
            
            %Get list of supported views, converted to corresponding classes
            vws = cellfun(@func2str, mgr.SupportedViews(:,2), 'UniformOutput', false);
            
            %Add to list
            cls(end+1:end+numel(vws)) = vws;
            
            %Get metaclass details whilst we're here
            mcs = cellfun(@metaclass, cls, 'UniformOutput', false);
            
            %Down-select to mvc classes only
            cls(~strncmp('mvc.', cls, numel('mvc.'))) = [];
            
            %Where are these classes ?
            cls = cellfun(@which, cls, 'UniformOutput', false);
            
            %But in the deployment folder, not the original
            cls = strrep(cls, src, dd);
            
            %Anything within the mvc. package is a candidate or removal
            fdd = fdd(strncmp(fullfile(dd, '+mvc'), fdd, numel(fullfile(dd, '+mvc'))));
            
            %But only if it is also NOT on this list of required classes
            fdd = fdd(~ismember(fdd, cls));
            
            %Remove them
            delete(fdd{:});
            
            %Strip references to Debuggable and Deployable from any classes in the deploy folder
            
        end

        function deployToCompilerProject(obj, dd, varargin)
        
            %Pass it on
            applicationCompiler(obj.EntryPointFcn);
            
        end

        function deployToMATLABApp(obj, dd, varargin)
        
            %Correspnding prj file ?
            [~,n,~] = fileparts(obj.EntryPointFcn);
            prj = [n, '.prj'];
            if exist(prj, 'file') == 2
                
                %Pass it on
                matlab.apputil.create(obj.EntryPointFcn);
        
            else
                
                %Create new
                matlab.apputil.create;
                
            end
            
        end

        function deployToExecutable(obj, dd, varargin)
        
            %Might take a while
            disp('Creating standlone executable...please wait');
            
            %Pass it on to compiler
            mcc('-e', obj.EntryPointFcn, '-d', dd);
            
        end
    
    end
    
end
