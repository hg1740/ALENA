classdef (ConstructOnLoad) Application ... % Provides generic functionality required of a "typical" Application
        < mvc.mixin.Serializable ...       % including ability to persist a session (i.e. save / load to file)
        & mvc.mixin.Importable ...         % and the ability to import content from some external file format
        & mvc.mixin.Exportable ...         % and the ability to export content to some external file format
        & mvc.mixin.Viewable ...           % and the ability to view content, and to manage multiple views
        & mvc.mixin.Documentable ...       % and the ability to provide documentation to the user
        & mvc.mixin.Deployable             % and the ability for the Application to be deployed in various ways
    
    methods % Construction
        
        function obj = Application(varargin)
            
            %Anything in ?
            if ~isempty(varargin)
                
                %Apply
                set(obj, varargin{:});
                
            end
            
            %If caller does not want it back
            if nargout == 0
                
                %Open in corresponding view manager
                mvc.view.ViewManager(obj);
                
            end
            
        end

    end
    
end
