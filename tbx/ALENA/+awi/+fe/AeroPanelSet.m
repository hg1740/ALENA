classdef AeroPanelSet < awi.fe.FEBaseClass
    %AeroPanelSet Defines a list of aerodynamic panels.
    %
    % The definition of the 'AeroPanelSet' object matches that of the
    % AELIST1 bulk data type from MSC.Nastran.
    
    %Primary properties
    properties
        %ID number of this aerodynamic set
        SID
        %ID number of the aerodynamic panels in this set
        E
    end
    
    %Store a reference to the 'awi.fe' objects
    properties (Hidden = true)
        %Handle to the 'awi.fe.AeroPanel' objects
        AeroPanels
    end
    
    methods % set / get
        function set.SID(obj, val)  %set.SID 
            %set.SID Set method for the property 'SID'.
            
            %Pass it on to the underlying 'ID' property.
            obj.ID = val;
            
        end
        function set.E(obj, val)    %set.E   
            %set.E Set method for the property 'E'.
            
            validateattributes(val, {'numeric'}, {'column', 'integer'}, ...
                class(obj), 'E');
            obj.E = val;
        end
        function val = get.SID(obj) %get.SID 
            %get.SID Get method for property 'SID'.
            val = obj.ID;
        end
        function val = get.E(obj)   %get.E   
            %get.E Get method for the property 'E'.
            %
            % If the object has been assigned a handle to its
            % 'awi.fe.AeroPanel' objects then always use their ID number,
            % else use 'E'.
            if isempty(obj.AeroPanels)
                val = obj.E;
            else
                temp = arrayfun(@(ap) (ap.EID : ap.EID + ap.NumPanels - 1), ...
                    obj.AeroPanels, 'Unif', false);
                val = horzcat(temp{:});
            end
        end
    end
    
    methods % construction
        function obj = AeroPanelSet
            
            %Make a note of the property names
            addFEProp(obj, 'SID', 'E');
            
        end
    end
    
    methods % writing data in MSC.Nastran format
        function writeToFile(obj, fid, bComment)
            %writeToFile Write the data for the 'awi.fe.AeroPanelSet'
            %object into a text file using the format of the MSC.Nastran
            %'AELIST' bulk data entry.
            
            %By default, do not close the file
            bClose = false;
                        
            if nargin < 2 %Ask the user for the file
                fName = awi.fe.FEBaseClass.getBulkDataFile;
                bClose = true;
                fid = fopen(fName, 'w');                
            end
            
            if nargin < 3 %Comments by standard
                bComment = true;
            end
                        
            if bComment %Helpful comments?
                comment = 'AELIST : Defines a list of aerodynamic elements.';
                awi.fe.FEBaseClass.writeComment(comment, fid);                  
            end
            
            awi.fe.FEBaseClass.writeColumnDelimiter(fid, '8');
                
            %How many objects?
            nObj = numel(obj);
                        
            %Plot in a loop for now...            
            for ii = 1 : nObj 
                %TODO - Vectorise this!
                awi.fe.FEBaseClass.writeListDataCard(fid, 'AELIST', obj(ii), {'E'});                
            end
            
            if bClose %Close the file?
                fclose(fid);
            end
            
        end
    end    
    
end

