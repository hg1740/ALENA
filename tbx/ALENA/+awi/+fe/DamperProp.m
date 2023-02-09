classdef DamperProp < awi.fe.ScalarProp
    %DamperProp Describes the damping coefficient a generic damper element.
    %
    % The defintion of the 'DamperProp' object matches that of the PDAMP
    % bulk data type from MSC.Nastran.
    
    %Primary properties
    properties
        %Damping coefficient (force per unit velocity)
        B
    end
    
    methods % set / get
        function set.B(obj, val)    %set.B   
            obj.PropVal = val;
        end       
        function val = get.B(obj)   %get.B
           val = obj.PropVal; 
        end
    end
    
    methods % construction
        function obj = DamperProp
            
            addFEProp(obj, 'PID', 'B');
            
        end
    end
    
    methods % writing data in MSC.Nastran format
        function writeToFile(obj, fid, bComment)
            %writeToFile Write the data for the 'awi.fe.DamperProp' object
            %into a text file using the format of the MSC.Nastran PDAMP
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
                comment = ['PDAMP : Specifies the damping value of a ', ...
                    'scalar damper element'];
                awi.fe.FEBaseClass.writeComment(comment, fid);                   
            end
            
            awi.fe.FEBaseClass.writeColumnDelimiter(fid, '16');
            
            nObj = numel(obj);
            
            %Card name
            nam   = repmat({'PDAMP*'}, [1, nObj]);
            
            %Set up the format for printing
            data = [nam ; {obj.ID} ; {obj.B}];
            
            %Write in 16-character column width as standard
            format = '%-8s%-16i%#-16.6g\r\n';
            
            %Write the data to the file
            fprintf(fid, format, data{:});
            
            if bClose %Close the file?
                fclose(fid);
            end
            
        end
    end
    
end

