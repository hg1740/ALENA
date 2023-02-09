classdef (ConstructOnLoad) Auditable < matlab.mixin.SetGet
    
    properties (SetAccess = protected, SetObservable)
        
        %TODO- Do we actually need to call "rebuild" whenever the
        %AuditTrail is modified? Is this prohibitively expensive?
        
        %Storage for audit trail
        AuditTrail;
        
        %What field to access from object, to treat as the name associated with audit trail entry
        AuditTrailField = 'Name';
        
    end
    
    properties (Dependent)
        
        AuditTrailExpanded;
        
    end
    
    properties (Transient, Access = private)
        
        %Placeholder for listeners
        AuditTrailListeners;
        
    end
    
    methods % get/set
        
        function val = get.AuditTrailExpanded(obj)
    
            %Start here
            val = obj.AuditTrail;
            
            %Any content in form of a cell ?
            b = cellfun(@iscell, val.What);
            
            %For each in turn (in reverse order)
            for i = flip(find(b(:).'))
               
                %Split into section above, at and below the line of interest
                valu = val(1:i-1,:);
                valb = val(i+1:end,:);
                
                %Expand line of interest
                valm = val(i,:);
                dat = valm.What{1};
                for j = 1:numel(dat)
                    valm.What{j, 1} = dat{j}; %Edit C.Szczyglowski - Need to force column vector so that vertcat works correctly
                end
                
                %Assign single "When", "User", "Object" to all - Edit
                %C.Szczyglowski - Do we need to repeat all values or is the
                %expansion implicit? The text data in the file looks much
                %neater if we replicate the data...
                valm.When   = repmat(valm.When(1)  , size(valm.What));
                valm.User   = repmat(valm.User(1)  , size(valm.What));
                valm.Object = repmat(valm.Object(1), size(valm.What));
                
                %Combine
                val = [valu; valm; valb]; %Edit C.Szczyglowski - Variable 'valb' was not being used.
                
            end
            
        end
        
    end
    
    methods % constructor/destructor
            
        function obj = Auditable
            
            %Initialise content
            obj.AuditTrail = obj.makeAuditTrailEntry;
            
            %Extend proprety group detail (if applicable)
            if isa(obj, 'mvc.mixin.Nameable')
                
                %Add as new property group
                addPropertyGroup(obj, 'Audit Trail', 'AuditTrail', 'Audit trail detail associated with this instance');
                
                %But do NOT want the audit trail to appear in exported content
                setPropertyAttribute(obj, 'Audit Trail', 'AuditTrail', 'Export', false);
                
            end
            
            %Extend context (if applicable)
            if isa(obj, 'mvc.mixin.Contextable')
                
                %Specify a low priority, so the "Audit" menu appears towards end
                addContext(obj, 10, 'Audit trail', '');
                
                %Always include options to show and clear
                addContext(obj, ...
                    'Audit trail>Show', 'showAuditTrail', ...
                    'Audit trail>Clear', 'clearAuditTrail');
                
                %For test / debug purposes (when released, may not want to make these too easily accessible)
                addContext(obj, ... 
                    'Audit trail>|Start', 'startAuditing', ...
                    'Audit trail>Stop', 'stopAuditing');
                
            end

        end
        
    end
    
    methods % manipulating the audit trail details
        
        function stopAuditing(obj, varargin)
        
            %If nothing specific passed in
            if nargin == 1
                
                %If this object is Collectable, and has content in its collection
                if isa(obj, 'mvc.mixin.Collectable') && obj.HasChildren
                    
                    %Applicable to whom ?
                    b = arrayfun(@(x)isa(x, 'mvc.mixin.Auditable'), obj.Children);
                    
                    %Pass it on
                    arrayfun(@stopAuditing, obj.Children(b));
                    
                end
                
            end
            
            %Anything to do ?
            if isempty(obj.AuditTrailListeners)
                
                %No
                return;
            
            elseif isempty(varargin)
                
                %Stop listening to absolutely everything
                delete(obj.AuditTrailListeners);
                obj.AuditTrailListeners = [];
                
            elseif iscellstr(varargin)
                
                %Stop listening to specified properties
                [b, idx] = ismember(arrayfun(@(x)x.Source{1}.Name, obj.AuditTrailListeners, 'UniformOutput', false), varargin);
                delete(obj.AuditTrailListeners(idx(b)));
                obj.AuditTrailListeners(idx(b)) = [];
                
            end            
            
        end
        
        function startAuditing(obj, varargin)
        
            %If nothing specific passed in
            if nargin == 1
                
                %If this object is Collectable, and has content in its collection
                if isa(obj, 'mvc.mixin.Collectable') && obj.HasChildren
                    
                    %Applicable to whom ?
                    b = arrayfun(@(x)isa(x, 'mvc.mixin.Auditable'), obj.Children);
                    
                    %Pass it on
                    arrayfun(@startAuditing, obj.Children(b));
                    
                end
                
                %Audit everything
                varargin = allProperties(obj);
                
                %But you wouldn't want to maintain an audit trail of the audit trail - that would be plain daft
                varargin = setdiff(varargin, 'AuditTrail');
                
            end
            
            %Properties listed in varargin are rendered auditable
            %TODO: figure out how to capture old and new value of property within audit trail
            % (can this be done with a single listener, or do we need to listen to preset AND postset ?)
            obj.AuditTrailListeners = addlisteners(obj, varargin, @i_addAuditTrailEntry);
            
            function i_addAuditTrailEntry(prp, varargin)
               
                %Pass it on
                addAuditTrailEntry(obj, [prp.Name, ' changed']);
                
            end
            
        end
        
        function addAuditTrailEntry(obj, varargin)
            
            %Create entry
            auditTrailEntry = obj.makeAuditTrailEntry(varargin{:});
            
            %Prepend to store
            obj.AuditTrail = [auditTrailEntry; obj.AuditTrail];

        end
                
        function clearAuditTrail(obj, bNoPrompt)
        
            %Caller can force no-prompt
            if nargin > 1 && bNoPrompt
                
            elseif ~confirm(obj) %TODO: May need to control access to this, e.g. with a 'Locked' flag ??
                
                %Bail out
                return;
                
            end
            
            %Easy
            obj.AuditTrail = obj.makeAuditTrailEntry;
            
        end
        
        function fn = showAuditTrail(obj, fn)
        
            %Where ?
            if nargin < 2
                
                %Make something up
                fn = [tempname, '.txt'];
                
            end
            
            %Send table to destination file
            writetable(obj.AuditTrailExpanded, fn, 'Delimiter', '\t');
            
            %But we need some basic header info too, I think, so get the content back again
            C = fileread(fn);
            
            %Prepend header info
            C = [{sprintf('%s\nAudit trail exported by user %s at %s\n\n', dlgtitle(obj), getenv('UserName'), datetime)}; C];
            
            %And write back again
            fid = fopen(fn, 'w');
            assert(fid ~= -1, 'failed to open file');
            fprintf(fid, '%s\n', C{:});
            fclose(fid);
            
            %If caller did not want the filename back
            if nargout == 0
                
                %Show us
                edit(fn);
                
            end
            
        end
        
    end
    
    methods (Access = private)
    
        function auditTrailEntry = makeAuditTrailEntry(obj, what, who, when)
            %
            % Helper function ensures consistency how-so-ever an entry in Audit Trail is made
             
            %Might have nothing to do
            if nargin > 1 && istable(what) && all(strcmp(what.Properties.VariableNames, {'When', 'User', 'Object', 'What'}))
                
                %Just send it straight back
                auditTrailEntry = what;
                
                %Caller may overwrite who and when
                if nargin > 2
                    auditTrailEntry.User = who;
                end
                if nargin > 3
                    auditTrailEntry.When = when;
                end
                
            else
                
                %Name for this entry ?
                name = obj.(obj.AuditTrailField);
                
                %Sensible defaults
                if nargin < 4
                    when = datetime;
                end
                if nargin < 3
                    who = getenv('Username');
                end
                
                %If no 'what' is provided
                if nargin < 2
                    
                    %We are just initialising an empty table
                    bEmpty = true;
                    what = [];
                    
                else
                    bEmpty = false;
                end
                
                %Create a one-line table
                auditTrailEntry = cell2table({when, who, name, what}, 'VariableNames', {'When', 'User', 'Object', 'What'});
                %
                % Or better still, a time-table
                %auditTrailEntry = timetable(when, {who}, {name}, {what}, 'VariableNames', {'User', 'Object', 'What'});
                
                %If initialising to empty
                if bEmpty
                    auditTrailEntry(1,:) = [];
                end
            
            end
            
        end
        
    end
    
end
