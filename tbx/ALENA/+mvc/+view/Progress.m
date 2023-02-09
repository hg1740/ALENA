classdef Progress < mvc.view.Container
%Progress view is just a scrolling listbox of lines of text,
% used to present progress during some generic sequence of processes

    properties
        
        %The model to which we are attached
        Model;
        
        %Function to call on selection change
        SelectionChangeFcn;
        
        %User may choose to cancel during whatever is being reported
        Cancelled = false;
        
    end
    
    properties (Dependent)
        
        %Access to the string in the list
        String;
        Value;
        
        %Track currently selected item in model
        %   - TODO : Remove this
        %Selection;
        
    end
    
    properties
       Selection 
    end
    
    properties (Access = protected)
        
        %The listbox in which we display the text
        hText;
        
        %A Cancel button
        hCancel;
        
    end
    
    methods
        
        function val = get.String(obj)
            
            %Pass it on
            val = get(obj.hText, 'String');
            
            %Ensure return as cell
            if isempty(val)
                val = {};
            elseif ischar(val)
                val = {val};
            else
                assert(iscell(val), 'bad data');
            end
            
        end
        
        function set.String(obj, val)
            
            %Ensure set as cell
            if isempty(val)
                val = {};
            elseif ischar(val)
                val = {val};
            else
                assert(iscell(val), 'invalid data');
            end
            
            %Pass it on
            set(obj.hText, 'String', val);
            
        end
        
        function val = get.Value(obj)
            
            %Pass it on
            val = get(obj.hText, 'Value');
            
        end
        
        function set.Value(obj, val)
            
            %Ensure valid ?
            
            %Pass it on
            set(obj.hText, 'Value', val);
            
        end
        
        function set.Selection(obj, val)
            
            %This view doesn't care about the selection,
            % so do nothing
            
        end      
        
        function obj = Progress(model, varargin )
            
            %To allow the constructor to be called with no arguments (e.g. for test purposes)
            if nargin == 0
                
                %Make something up
                model = mvc.model.DrawableThing.defaultTree;
                
            end
            
            %Call superclass constructor
            obj@mvc.view.Container( varargin{:} );
            
            %Create a VBox
            hv = uiextras.VBox('Parent', obj.UIContainer, 'Spacing', 6, 'Padding', 6);
            
            %Add a listbox to contain the text
            obj.hText = uicontrol('Parent', hv, ...
                'Units', 'normalized', ...
                'Position', [0, 0, 1, 1], ...
                'Style', 'listbox', ...
                'Min', 0, ...
                'Max', 2);
            
            %Add a Cancel button
            hb = uiextras.HButtonBox('Parent', hv);
            hv.Heights(end) = 50;
            obj.hCancel = uicontrol('Parent', hb, ...
                'String', 'Cancel', ...
                'Callback', @obj.cancel);
            
            %Store the model
            obj.Model = model;
            
        end
        
        function [bOK, strs] = update(obj, str)
            
            %Return flag indicating existence of listbox, and state of Cancel button
            bOK = ishghandle(obj.hText) && ~obj.Cancelled;
            
            %Not ideal, but need a way to allow caller to explicitly close
            bClose = false;
            
            %Anything to add ?
            if ishghandle(obj.hText) && nargin > 1 && ~isempty(str)
                
                %Allow explicit close
                if strcmpi(str, 'close')
                    
                    %Make a note for later
                    bClose = true;
                    
                elseif isstrprop(str(1), 'punct')
                    
                    %A little presumptious, perhaps, but if the string starts with a punctuation character
                    % append existing last-line of listbox
                    obj.String{end} = [obj.String{end}, str];
                    
                else
                    
                    %Add new line item to list
                    obj.String{end+1} = str;
                    obj.Value = numel(obj.String);
                    
                end
                
                %Force an update
                drawnow;
                
            end
            
            %Caller might want the entire message back too
            if nargout > 1
                
                %Pass it on
                strs = obj.String;
                                
            end
            
            %But need some way of allowing caller to close - this is not ideal
            if bClose
                
                %Try to detach the 'docked' view if it is parented by a
                %view manager
                try
                    %Find the figure
                    hFig = ancestor(obj, 'Figure');
                    if isempty(hFig) %Escape route
                        obj.delete
                    end
                    %Check for 'UserData' and extract the
                    %'mvc.view.ViewManager' object.
                    ud = hFig.UserData;
                    if iscell(ud)
                        idx = cellfun(@(x) isa(x, 'mvc.view.ViewManager'), ud);
                        vwm = ud{idx};
                    else
                        idx = arrayfun(@(x) isa(x, 'mvc.view.ViewManager'), ud);
                        vwm = ud(idx);
                    end
                    if isempty(vwm) %Escape route
                        obj.delete
                    end
                    %Detach the view then delete it
                    obj = vwm.removeView(obj.Title);
                    delete(obj);
                catch                    
                    %Just delete it
                    obj.delete;
                end
            end
            
        end
        
    end
    
    methods (Access = protected)
        
        function cancel(obj, varargin)
            
            %Make a note
            obj.Cancelled = true;
            
        end

    end
    
end
