classdef (ConstructOnLoad) Simulatable < matlab.mixin.SetGet
    %
    % Simulatable provides functionality associated with Simulink
    
    properties (AbortSet, SetObservable)
    
        %The model itself
        Model; % Hmmm, any danger of confusion over terminology here ?
        
        %To configure the run
        SignalLoggingName = 'logsout';
        
    end    
    
    properties (Dependent)
        
        ModelHandle;
        ModelFileName;
        
    end
    
    methods % get/set
        
        function set.ModelFileName(obj, val)
            
            %Split
            [p,n,e] = fileparts(val);
            
            %Anything ?
            if isempty(p)
                
                %Get whatever
                [p, n, e] = fileparts(which([n, e]));
                assert(~isempty(p), 'file not found');
                
            end
            
            %Open the model
            open_system([n,e]);
            sys = bdroot;
            
            %Make sure we're good to go
            assert(strcmpi(fullfile(p,[n,e]), get_param(sys, 'FileName')), 'path mismatch');
            
            %Make a note
            obj.Model = get_param(sys, 'Name');            
            
        end
        
        function val = get.ModelHandle(obj)
    
            %Anything yet ?
            if isempty(obj.Model)
                
                %No
                val = [];
                
            else
                
                %Look for it
                val = find_system('FindAll', 'on', ... % So we return a handle
                    'FirstResultOnly', 'on', ... % because there can be only one
                    'Type', 'block_diagram', ...
                    'Name', obj.Model);
                            
            end
            
        end
        
        function val = get.ModelFileName(obj)
    
            %Start from handle
            h = obj.ModelHandle;
            
            %Anything yet ?
            if isempty(h)
                
                %No
                val = [];
                
            else
                
                %Ask it
                val = get_param(h, 'FileName');
                
            end
            
        end
        
    end
    
    methods % construction / destruction
        
        function obj = Simulatable(varargin)
            
            %Extend property groups (if applicable)
            if isa(obj, 'mvc.mixin.Nameable')
                obj.addPropertyGroup('Simulink', ...
                    'Model', 'Model', ...
                    'SignalLoggingName', 'Signal logging name');
            
                %Tweak naming convention
                obj.DefaultField = 'Model';
            
            end
            
            %Extend context (if applicable)
            if isa(obj, 'mvc.mixin.Contextable')
                addContext(obj, ...
                    'Simulink>Model...', 'onModel', ...
                    'Simulink>Open...', 'onOpen', ...
                    'Simulink>Run...', 'onRun', ...
                    'Simulink>Sweep...', 'onSweep');
 
            end
            
        end
        
    end
    
    methods % thin wrappers around standard Simulink
    
        function varargout = find_system(obj, varargin)
            
            %Pass it on
            [varargout{1:nargout}] = find_system(obj.ModelHandle, varargin{:});
            
        end
        
        function varargout = get_param(obj, varargin)
            
            %Pass it on
            [varargout{1:nargout}] = get_param(obj.ModelHandle, varargin{:});
            
        end
        
        function varargout = set_param(obj, varargin)
            
            %Pass it on
            [varargout{1:nargout}] = set_param(obj.ModelHandle, varargin{:});
            
        end
        
        function varargout = open_system(obj, varargin)
            
            %Pass it on
            [varargout{1:nargout}] = open_system(obj.ModelHandle, varargin{:});
            
        end
        
        function res = sim(obj, varargin)
            
            %Additional args
            if ~isempty(obj.SignalLoggingName)
                args = {'SignalLogging', 'on', 'SignalLoggingName', obj.SignalLoggingName};
            else
                args = {};
            end
            
            %Pass it on
            res = sim(obj.Model, args{:}, varargin{:});
            
            %If caller did not want the result back
            if nargout == 0
                
                %Add to internal store
                obj.add(mvc.model.SimulationResult('Result', res));
                
            end
            
        end
        
    end
    
    methods % interaction with the associated Simulink model
       
        function onSweep(obj, varargin)
            
            %Sweep what ?
            val = inputdlg(obj, {5}, {'How many runs ?'});
            if isempty(val)
                return;
            end
            
            %Unpick
            n = val{1};
            
            %Create new sweepset
            swp = obj.add(mvc.model.SimulationSweep);
            
            %For each member of sweep
            for i = 1:n
                
                %Run the model
                res = sim(obj, varargin{:});
                
                %Add result to sweepset
                swp.add(mvc.model.SimulationResult('Result', res));                
                
            end
            
        end
        
        function onRun(obj, varargin)
            
            %Anything yet ?
            if isempty(obj.Model)
                if isempty(obj.onOpen)
                    return;
                end
            end
            
            %Pass it on
            sim(obj);
            
        end
       
        function sys = onOpen(obj, varargin)
            
            %Anything yet ?
            if isempty(obj.Model)
                
                %Pass it on
                sys = obj.onModel;
                if isempty(sys)
                    return;
                end
                
            end
            
            %Pass it on
            open_system(obj);
            sys = bdroot;
            
        end
        
        function sys = onModel(obj, varargin)
        
            %Ask the user
            [f,p] = uigetfile({'*.slx', 'Simulink model files (*.slx)'}, 'Choose model...', char(obj.Model));

            %Cancelled ?
            if isempty(f) || isequal(f, 0)
                sys = [];
                return;
            end
            
            %Pass it on
            obj.ModelFileName = fullfile(p,f);
            
        end
        
        
    end    
    
end
