classdef (ConstructOnLoad) Analysable < matlab.mixin.SetGet
    %
    % The Analysable class provides functionality that allows an object
    % to know about, and access, applicable analyses
    
    properties (SetAccess = protected, Transient, AbortSet)
        
        %What analyses do we support ?
        SupportedAnalyses = struct('Name', {}, ...
            'Description', {}, ...
            'EntryPointFcn', {}, ...
            'ApplicableFcn', {}, ...
            'AdditionalArgs', {});
        
    end
    
    properties (Dependent)
    
        %To help exploit contextable functionality
        SupportedAnalysesContextStruct;
        
    end
    
    methods % get/set
        
        function val = get.SupportedAnalysesContextStruct(obj)
            
            %Create equivalent structure
            val = struct([]);
            for i = 1:numel(obj.SupportedAnalyses)
                
                %What's it called ?
                nam = obj.SupportedAnalyses(i).Name;
                
                %Is it enabled ?
                bEn = isempty(obj.SupportedAnalyses(i).ApplicableFcn) || obj.SupportedAnalyses(i).ApplicableFcn(obj);
                
                %Hence assign callback
                if bEn
                    cb = {obj.SupportedAnalyses(i).EntryPointFcn, obj.SupportedAnalyses(i).AdditionalArgs{:}};
                else
                    cb = [];
                end
                
                %Get help
                val = mvc.mixin.Contextable.addToContextStruct(val, nam, cb);
                
            end
            
        end
        
    end
    
    methods % construction / destruction
        
        function obj = Analysable(varargin)
            
            %Extend context (if applicable)
            if isa(obj, 'mvc.mixin.Contextable')
                addContext(obj, 'Analyse', 'analysemenu', ...
                    'Analyse>[not available]', []);
            end
                        
        end
        
    end
        
    methods
        
        function addAnalysis(obj, efcn, afcn, name, desc, varargin)
        
            %What have we got ?
            if nargin < 2 || isempty(efcn)
                
                %TODO: browse for file ??
                
            elseif isa(efcn, 'function_handle')
                
                %OK
                
            elseif ischar(efcn)
                
                %Cast to function handle
                efcn = str2func(efcn);
                
            else
                error('bad input');
            end
            
            %Applicability function - optional
            if nargin < 3
                afcn = [];
            end
            
            %Name ?
            if nargin < 4 || isempty(name)
                
                %Good guess
                name = func2str(efcn);
                
            end
            
            %Description ?
            if nargin < 5 || isempty(desc)
                
                %Good guess
                desc = help(which(func2str(efcn)));
                
            end

            %Build structure
            S = struct('Name', name, ...
                'Description', desc, ...
                'EntryPointFcn', efcn, ...
                'ApplicableFcn', afcn, ...
                'AdditionalArgs', {varargin});
            
            %Add to store
            obj.SupportedAnalyses(end+1) = S;
            
        end
        
    end
    
    methods % context menu maintenance
        
        function varargout = analysemenu(obj, varargin)
            
            %Requires help from contextable
            assert(isa(obj, 'mvc.mixin.Contextable'), 'MUST be contextable');
            
            %Need the calling handle (is this always safe ?)
            hc = gcbo;
            
            %Anything to worry about ?
            if isempty(hc)
                
                %No
                return;
                
            end
                
                %                 %Create equivalent structure
                %                 ctx = [];
                %                 for i = 1:numel(obj.SupportedAnalyses)
                %
                %                     %What's it called ?
                %                     nam = obj.SupportedAnalyses(i).Name;
                %
                %                     %Is it enabled ?
                %                     bEn = isempty(obj.SupportedAnalyses(i).ApplicableFcn) || obj.SupportedAnalyses(i).ApplicableFcn(obj);
                %
                %                     %Hence assign callback
                %                     if bEn
                %                         cb = {obj.SupportedAnalyses(i).EntryPointFcn, obj.SupportedAnalyses(i).AdditionalArgs{:}};
                %                     else
                %                         cb = [];
                %                     end
                %
                %                     %Get help
                %                     ctx = mvc.mixin.Contextable.addToContextStruct(ctx, nam, cb);
                %
                %                 end
            
                %Pass it on
                context(obj, hc, obj.SupportedAnalysesContextStruct);
                
                %             else
                %
                %                 %Need to cater separately for HG and Java objects
                %                 if ishghandle(hc)
                %
                %                     %Simple
                %                     hm = hc;
                %
                %                     %Clean sheet
                %                     delete(allchild(hm));
                %
                %                     %Get options
                %                     S = obj.SupportedAnalyses;
                %
                %                     %Anything ?
                %                     if isempty(S)
                %
                %                         %No
                %                         uimenu(hm, 'Label', '[none]', 'Enable', 'off');
                %
                %                     else
                %
                %                         %For each analysis
                %                         for i = 1:numel(S)
                %
                %                             %Is this function enabled ?
                %                             bEn = isempty(S(i).ApplicableFcn) || S(i).ApplicableFcn(obj);
                %
                %                             %Add to menu
                %                             uimenu(hm, 'Label', S(i).Name, ...
                %                                 'Enable', mvc.mixin.UiTools.bool2offon(bEn), ...
                %                                 'Callback', {@i_invoke, S(i)});
                %
                %                         end
                %
                %                     end
                %
                %                 elseif isa(hc, 'javahandle_withcallbacks.javax.swing.JMenu')
                %
                %                     %Clean sheet
                %                     jmenu = hc; % .getItem(idx);
                %                     for i = jmenu.getItemCount()-1:-1:0
                %                         jmenu.remove(i);
                %                     end
                %
                %                     %Get options
                %                     S = obj.SupportedAnalyses;
                %
                %                     %Anything ?
                %                     if isempty(S)
                %
                %                         %No - create a menu item with no action
                %                         jitem = javax.swing.JMenuItem('[none]');
                %                         jmenu.add(jitem);
                %
                %                     else
                %
                %                         %For each analysis
                %                         for i = 1:numel(S)
                %
                %                             %Create a menu item with action
                %                             jitem = javax.swing.JMenuItem(S(i).Name);
                %
                %                             %Is this function enabled ?
                %                             bEn = isempty(S(i).ApplicableFcn) || S(i).ApplicableFcn(obj);
                %                             if bEn
                %
                %                                 %Assign action
                %                                 set(jitem, 'ActionPerformedCallback', {@i_invoke, S(i)});
                %
                %                             else
                %
                %                                 %Disable
                %                                 jitem.setEnabled(false);
                %
                %                             end
                %
                %                             %And add to menu
                %                             jmenu.add(jitem);
                %
                %                         end
                %
                %                     end
                %
                %                 else
                %                     warning(['unhandled class ''', class(hc), '''']);
                %                 end
                %
                %             end
                %
                %             function i_invoke(hc, ~, S)
                %
                %                 %Might take some time
                %                 clu = mvc.mixin.UiTools.pointer(hc); %#ok<NASGU>
                %
                %                 %Careful
                %                 try
                %
                %                     %Pass it on
                %                     S.EntryPointFcn(obj, S.AdditionalArgs{:});
                %
                %                 catch err
                %
                %                     %What went wrong ?
                %                     uiwait(errordlg(obj, err, 'modal'));
                %
                %                 end
                %
                %             end
            
        end
        
    end
    
    methods (Sealed)
    
        function vw = analyse(obj, sel, varargin)
        
            %What are the options ?
            vws = obj.SupportedAnalyses;
            
            %Safety net
            assert(~isempty(vws), ['no supported analyses(s) for class ', class(obj)]);
            
            %View what ?
            if nargin < 2 || isempty(sel)
                
                %Any ambiguity ?
                if numel(vws) == 1
                    
                    %No
                    sel = 1;
                    
                else
                    
                    %Ask the user
                    sel = listdlg(obj, {vws.Name}, 'Choose analysis...');
                    
                    %Cancelled ?
                    if isempty(sel)
                        return;
                    end
            
                end
                
            elseif isnumeric(sel)
                
                %Get by index
            
            elseif ischar(sel) || iscellstr(sel)
                
                %Look for match
                %For backward compat with 2015b, do not use "contains"
                % sel = find(contains({vws.Name}, sel, 'IgnoreCase', true));
                sel = find(~cellfun(@(x)strfind(lower(sel), lower(x)), {vws.Name}));
                assert(~isempty(sel), 'no match');
                assert(numel(sel) == 1, 'ambiguous match');
                
            else
                error('bad input');
            end
            
            %Allow for multiple selection
            for i = 1:numel(sel)
                vw(i) = vws(sel(i)).EntryPointFcn(obj, vws(sel(i)).AdditionalArgs, varargin{:});
            end
            
        end
        
    end
    
end
