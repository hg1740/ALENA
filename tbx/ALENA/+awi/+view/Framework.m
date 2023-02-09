classdef Framework < mvc.view.ViewManager
    
    methods
        
        function obj = Framework(model, varargin)
            
            %To allow the constructor to be called with no arguments (e.g. for test purposes)
            if nargin == 0
                
                %Make something up
                model = awi.model.Framework.defaultFramework;
                
            end
                
            %Pass it on
            obj@mvc.view.ViewManager(model, varargin{:});
            
            % Extend list of supported views
            obj.SupportedViews(end+1,:) = {'Static Results View', @awi.view.BeamResultsViewer};
            obj.SupportedViews(end+1,:) = {'Mass Distribution'  , @awi.view.MassDistribution};
            
            %Good name
            obj.GuiName = model.Name;
            
            %If the arrangement of views has NOT been applied by base class
            if isempty(obj.Views)
                
                %Start with a Tree View on the left-hand panel
                obj.addView('Tree');
                
                %Properties and Drawing view on right
                obj.addPanel;
                obj.PanelWidths = [-1 -4];
                obj.addView('Properties');
                obj.addView('Drawing');
            
            end
            
        end
        
    end
    
end

