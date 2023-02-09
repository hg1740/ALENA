classdef (ConstructOnLoad) Viewable < matlab.mixin.SetGet
    %
    % The Viewable class provides functionality that allows an object
    % to know about, and access, suitable view classes
    
    properties (SetAccess = protected, Transient, AbortSet)
        
        %What view(s) do we support ?
        SupportedViews = {};
        
    end
    
    events ( NotifyAccess = private )
        ViewContext;                      
    end
    
    properties (Dependent)
        
        DefaultView;
        
    end
    
    methods % get/set
    
        function val = get.DefaultView(obj)
                       
            %Start with class name, tokenised into package/subpackage etc
            tok = strsplit(class(obj), '.');
            
            %If only one token
            if numel(tok) == 1
                
                %There's nothing much we can do
                val = [];
                
            else
                
                %Replace the penultimate token with the name 'view'
                tok{end-1} = 'view';
                
                %Rebuild
                val = strjoin(tok, '.');
                
                %Does this class exist ?
                if exist(val, 'class') == 8
                    
                    %Yes - convert to function handle
                    val = str2func(val);
                    
                else
                    
                    %No - can't help
                    val = [];
                    
                end
                
            end
            
        end
    
        function val = get.SupportedViews(obj)
                  
            %Start here
            val = obj.SupportedViews;
            
            %If default view exists, and is NOT on the list of supported views
            def = obj.DefaultView;
            
            %Combine, ensuring that the default view always appears at top of list
            val = [{def}; val];
            
            %But cater for possibility that user has manually added the default too
            [~, udx] = unique(cellfun(@func2str, val, 'UniformOutput', false));
            val = val(udx);
                        
        end
        
    end
    
    methods % construction / destruction
        
        function obj = Viewable(varargin)
            
            %Extend context (if applicable)
            if isa(obj, 'mvc.mixin.Contextable')
                addContext(obj, ...
                    'View>New...', 'view', ...
                    'View>Close...', 'viewClose');                
                
                %Slightly experimental
                addContext(obj, ...
                    'View', 'viewContext', ...
                    'View>Add', []);
                
            end
                        
        end
        
        function varargout = viewContext(obj, varargin)
        
            %We don't actually do much here - just raise the event
            % to which others may choose to respond
            notify(obj, 'ViewContext');
            
            disp('view context...');
            
            %Find a View manager associated with this object
            
            %Pass control of contexts on
            
        end
        
    end
    
    methods (Sealed)
    
        function vw = view(obj, sel, varargin)
        
            %What are the options ?
            vws = obj.SupportedViews;
            
            %Safety net
            assert(~isempty(vws), ['no supported view(s) for class ', class(obj)]);
            
            %View what ?
            if nargin < 2 || isempty(sel)
                
                %Any ambiguity ?
                if numel(vws) == 1
                    
                    %No
                    sel = 1;
                    
                else
                    
                    %Ask the user
                    sel = listdlg(obj, cellfun(@func2str, wvs, 'UniformOutput', false), ...
                        'Choose view...');
                    
                    %Cancelled ?
                    if isempty(sel)
                        return;
                    end
            
                end
                
            elseif isnumeric(sel)
                
                %Get by index
            
            elseif ischar(sel) || iscellstr(sel)
                
                %Look for match
                [b, sel] = ismember(sel, cellfun(@func2str, wvs, 'UniformOutput', false));
                assert(all(b), 'no match');
                
            else
                error('bad input');
            end
            
            %Allow for multiple selection
            for i = 1:numel(sel)
                vw(i) = vws{sel(i)}(obj, varargin{:});
            end
            
        end
        
    end
    
end
