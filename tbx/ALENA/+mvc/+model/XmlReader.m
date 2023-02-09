classdef (ConstructOnLoad) XmlReader < mvc.model.Application ... is an Application
        & mvc.model.XmlItem % and it is an XML item (i.e. the root node of an XML doc)
    
    methods % Construction
        
        function obj = XmlReader(varargin)
            
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
