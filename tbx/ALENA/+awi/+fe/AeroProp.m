classdef AeroProp < awi.fe.FEBaseClass
    %AeroProp Defines associated bodies for the panels in the 
    %Doublet-Lattice method.
    %
    % The definition of the 'AeroProp' object matches that of the PAERO1
    % bulk data type from MSC.Nastran.
    
    %Primary Properties
    properties
        %Identification number for this aerodynamic property. 
        PID
        %Identification number of associated aerodynamic bodies
        B
    end
    
    properties
        %Handle to the ...
        AeroBodies
    end
    
    methods % set / get
        function set.PID(obj, val)  %set.PID 
            %set.PID Set method for the property 'PID'.
            %
            % Passes the value straight to the inherited 'ID' property
            % which will validate and store the ID number.
            obj.ID = val;
        end
        function val = get.PID(obj) %get.PID 
            %get.PID Get method for the property 'PID'.
            val = obj.ID;
        end
        function val = get.B(obj) %get.B
            %get.B Get method for the property 'B'.
            if isempty(obj.AeroBodies)
                val = obj.B;
            else
                val = [obj.AeroBodies.ID];
            end
        end
    end
    
    methods % construction
        function obj = AeroProp
            
            %Make a note of the property names
            addFEProp(obj, 'PID', 'B');
            
        end
    end
    
    methods % writing data in MSC.Nastran format
        function writeToFile(obj, fid, bComment)
            %writeToFile Write the data for the 'awi.fe.AeroProp' object
            %into a text file using the format of the MSC.Nastran 'PAERO1'
            %bulk data entry.           
            
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
                comment = ['PAERO1 : Defines associated bodies for ', ...
                    'the panels in the Doublet-Lattice method.'];
                awi.fe.FEBaseClass.writeComment(comment, fid);                   
            end
            
            awi.fe.FEBaseClass.writeColumnDelimiter(fid, '8');
               
            %How many objects?
            nObj = numel(obj);
                        
            %Plot in a loop for now...            
            for ii = 1 : nObj 
                %TODO - Vectorise this!
                awi.fe.FEBaseClass.writeListDataCard(fid, 'PAERO1', obj(ii), {'B'});                
            end
            
            if bClose %Close the file?
                fclose(fid);
            end
            
        end
    end
    
end

