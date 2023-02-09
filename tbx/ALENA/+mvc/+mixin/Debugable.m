classdef (ConstructOnLoad) Debugable < matlab.mixin.SetGet
    %
    % Debugable provides debug functionality useful to developers:
    %
    %  obj2base: send a copy of the current object into the base workspace.
    %  codeedit: open source code associated with current object in MATLAB editor.
    
    methods % construction / destruction
        
        function obj = Debugable(varargin)
            
            %Extend context (if applicable)
            if isa(obj, 'mvc.mixin.Contextable') && ~isdeployed
                addContext(obj, Inf, ... % Specify a super-low priority, so this appears at end
                    '|Debug>Send to base...', 'obj2base', ...
                    '|Debug>Source...', 'codeedit', ...
                    '|Debug>Help...', 'codehelp', ...
                    '|Debug>Explore...', 'codefolder', ...
                    '|Debug>Profile...', 'profile', ...
                    '|Debug>DB stop...', 'dbstop', ...
                    '|Debug>Start over...', 'startAgain' ...
                    );
            end
            
        end
        
    end
    
    methods (Sealed)
        
        function startAgain(obj)
        
            %How do we start again ?
            if isa(obj, 'mvc.mixin.Deployable') && ~isempty(obj.EntryPointFcn)
                
                %Go with the defined entry point
                fcn = obj.EntryPointFcn;
                
            else
                
                %Just assume that we create a new instance of the applciation class
                fcn = class(obj);
                
            end
            
            %So the way we start again is with this
            str = ['close all force, clear all, clear classes, clc, ', fcn];
            
            %Do it
            evalin('base', str);
            
        end
        
        function dbstop(obj, val)
            
            %Maintain a single (global !) state
            persistent DBSTOP;
            if isempty(DBSTOP)
                DBSTOP = 'off';
            end
            
            %Do what ?  Caller may pass it in
            if nargin < 2 || isempty(val)
                
                %But didn't do so.  What are the options ?
                options = {'on', 'off'};
                
                %Sensible default ?
                if strcmp(DBSTOP, 'on')
                    def = 'off';
                else
                    def = 'on';
                end
                
                %Ask the user
                val = choose(obj, 'DB stop if caught error ?', options{:}, def);
                
                %Cancelled ?
                if isempty(val)
                    
                    %Bail out
                    return;
                    
                end
            
            end
            
            %So do what
            if strcmpi(val, 'on')
                
                %Turn it on
                evalin('base', 'dbstop if caught error');
                
            elseif strcmpi(val, 'off')
                
                %Turn it off
                evalin('base', 'dbclear all');
                
            else
                %For time being
                error('invalid option');
            end
            
        end
        
        function profile(obj, val)
            
            %Maintain a single (global !) state
            persistent PROFILE;
            if isempty(PROFILE)
                PROFILE = 'off';
            end
            
            %Do what ?  Caller may pass it in
            if nargin < 2 || isempty(val)
                
                %But didn't do so.  What are the options ?
                options = {'on', 'off', 'report'};
                
                %Sensible default ?
                if strcmp(PROFILE, 'on')
                    def = 'report';
                else
                    def = 'on';
                end
                
                %Ask the user
                val = choose(obj, 'Profile ?', options{:}, def);
                
                %Cancelled ?
                if isempty(val)
                    
                    %Bail out
                    return;
                    
                end
            
            end
            
            %So do what
            profile(val);
                
        end
        
        function obj2base(obj, nam)
            
            %Using what variable name ?
            if nargin < 2 || isempty(nam)
                
                %Ask the user
                nam = inputdlg(obj, 'obj', ['Send instance ''', dlgtitle(obj), ''' to base workspace using what name ?']);
                
                %Cancelled ?
                if isempty(nam)
                    
                    %Bail out
                    return;
                    
                end

            end
            
            %Safety net
            nam = matlab.lang.makeValidName(nam);
            
            %Do it
            assignin('base', nam, obj);
            
        end
            
        function codehelp(obj, varargin)
           
            %Get list of possibilities
            cls = superclassList(obj);
            
            %Offer to user
            sel = listdlg(obj, cls);
            
            %Cancelled ?
            if isempty(sel)
                return;
            end
            
            %Down-select
            cls = cls(sel);
            
            %Go to doc centre, for each class
            cellfun(@doc, cls);
            
        end
            
        function codeedit(obj, varargin)
           
            %Get list of possibilities
            cls = superclassList(obj);
            
            %Offer to user
            sel = listdlg(obj, cls);
            
            %Cancelled ?
            if isempty(sel)
                return;
            end
            
            %Down-select
            cls = cls(sel);
            
            %Convert classname to code
            fn = cellfun(@which, cls, 'UniformOutput', false);
            
            %Edit it
            cellfun(@edit, fn);
            
        end
            
        function codefolder(obj, varargin)
           
            %Get list of possibilities
            cls = superclassList(obj);
            
            %Offer to user
            sel = listdlg(obj, cls);
            
            %Cancelled ?
            if isempty(sel)
                return;
            end
            
            %Down-select
            cls = cls(sel);
            
            %Convert classname to code
            fn = cellfun(@which, cls, 'UniformOutput', false);
            
            %Convert filename to folder
            dn = unique(cellfun(@fileparts, fn, 'UniformOutput', false));
            
            %Show us
            cellfun(@winopen, dn);
            
        end
        
        function varargout = superclassList(obj)
            %
            % Return a list of classes associated with this object and its super-classes
            
            %Start with metaclass for this object
            mc = metaclass(obj);
            
            %Pass it on to helper function, called recursively
            [varargout{1:nargout}] = i_classlist(mc);
            
            %Where the , looks like this
            function cls = i_classlist(mc, cls)
                
                %Anything yet ?
                if nargin < 2
                    
                    %No
                    cls = {};
                    
                end
                
                %Add details of this object
                cls(end+1:end+numel(mc)) = {mc.Name};
                
                %Superclasses
                for i = 1:numel(mc)
                    cls = i_classlist(mc(i).SuperclassList, cls);
                end
                
                %Not interested in the following
                cls(ismember(cls, {'handle', 'hgsetget'})) = [];
                cls(~cellfun(@isempty, strfind(cls, 'matlab.mixin.'))) = [];
                
            end
            
        end
        
    end
    
end
